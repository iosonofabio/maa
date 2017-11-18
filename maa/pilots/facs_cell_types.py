# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/17
content:    Try to see where in the sorting plots are successful and failed
            cells for different colon cell types (after RNA-Seq annotation).
'''
# Modules
import os
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset
from singlet.counts_table import CountsTable
from singlet.samplesheet import SampleSheet


# Globals
with open('facs_config.yml', 'r') as f:
    config = yaml.load(f)
    channels_to_genes = config['channels_to_genes']
    config = config['tissues']


# Functions
def parse_plate_metadata():
    import glob

    # Cache for faster access
    fn_cache = '../../data/plate_metadata/cache.tsv'
    if os.path.isfile(fn_cache):
        return pd.read_csv(
                fn_cache,
                sep='\t',
                index_col=0)

    fn_384 = glob.glob('../../data/plate_metadata/*384*.tsv')[0]
    md_384 = pd.read_csv(fn_384, sep='\t', index_col=0)
    md_384.index.name = 'name'
    md_384.rename(columns={
        'Experiment ID ': 'Experiment ID',
        },
        inplace=True)
    md_384['n.wells'] = 384

    fn_96 = glob.glob('../../data/plate_metadata/*96*.tsv')[0]
    md_96 = pd.read_csv(fn_96, sep='\t', index_col=0)
    md_96.index.name = 'name'
    md_96.rename(columns={
        'Mouse ID (age_#_sex)': 'mouse.id',
        'data.sorted': 'date.sorted',
        'EXP': 'Experiment ID',
        },
        inplace=True)
    md_96['n.wells'] = 96

    columns = [
            'date.sorted',
            'tissue',
            'subtissue',
            'mouse.id',
            'FACS.selection',
            'nozzle.size',
            'FACS.instument',
            'Experiment ID',
            'n.wells']

    md = pd.concat([md_384[columns], md_96[columns]], axis=0)

    # Kill non-plates
    md = md.reset_index().dropna().set_index('name')

    # Normalize subtissue
    sts = []
    for _, x in md.iterrows():
        st = x['subtissue']
        if (not st) or (str(st).lower() == 'nan') or (st.strip(' ?') == ''):
            sts.append('')
        else:
            sts.append(st.strip(' ').lower())
    md['subtissue'] = sts

    # Select only age 3
    age = []
    for _, x in md.iterrows():
        mid = x['mouse.id']
        if (not mid) or (str(mid).lower() == 'nan') or (mid.strip(' ?') == ''):
            age.append(-1)
        #NOTE: some diaphragm nonsense
        elif mid.startswith('HY_IP'):
            age.append(-1)
        # Some heart plates have more than one mouse (??)
        elif '&' in mid:
            age.append(-1)
        else:
            age.append(int(mid.split('_')[0]))
    md['mouse.age'] = age
    md = md.loc[md['mouse.age'] == 3]

    # Write cache
    md.to_csv(fn_cache, index=True, sep='\t')

    return md


def parse_facs_plate(tissue, plate):
    sorter = config[tissue].get('sorter', 'Sony')
    glb = config[tissue].get('facs_glob', '*')
    if sorter == 'Sony':
        return parse_facs_plate_sony(plate, glb=glb)
    else:
        return parse_facs_plate_aria(plate)


def parse_facs_plate_aria(plate, glb='*'):
    import glob
    import fcsparser

    out = {}

    fdn = '../../data/MACAFACS/index_fcs_{:}/'.format(glb)
    glob_fcs = fdn+plate+'.fcs'
    fn_fcs = glob.glob(glob_fcs)

    if len(fn_fcs) == 0:
        raise IOError('FCS file not found')
    if len(fn_fcs) > 1:
        raise IOError('Multiple FCS files found')
    fn_fcs = fn_fcs[0]

    glob_index = fdn+plate+'_Index.fcs'
    fn_index = glob.glob(glob_index)
    if len(fn_index) == 0:
        raise IOError('Index file not found')
    if len(fn_index) > 1:
        raise IOError('Multiple index files found')
    fn_index = fn_index[0]

    meta, data = fcsparser.parse(fn_fcs, reformat_meta=True)
    out['fcs_meta'] = meta
    out['fcs_data'] = data

    meta_index, data_index = fcsparser.parse(
            fn_index,
            meta_data_only=False,
            reformat_meta=True)

    # Figure in what wells the cells got sorted (the Aria's FCS is a mess)
    data_index['Index'] = 'A0'
    i = 1
    slstring = ''
    while 'INDEX SORTING LOCATIONS_'+str(i) in meta_index:
        slstring += meta_index['INDEX SORTING LOCATIONS_'+str(i)]
        i += 1

    itot = 0
    for sl in slstring.rstrip(';').split(';'):
        row, col = tuple(map(int, sl.split(',')))
        a24 = chr(65 + row)+str(col+1)
        data_index.loc[itot, 'Index'] = a24
        data_index.loc[itot, 'name'] = plate+'_'+a24
        itot += 1

    data_index.set_index('name', inplace=True)
    out['index_data'] = data_index
    return out


def parse_facs_plate_sony(plate, glb='*'):
    import glob
    import fcsparser

    out = {}

    fdn = '../../data/MACAFACS/index_fcs_{:}/'.format(glb)
    glob_fcs = fdn+'*'+plate+'*.fcs'
    fn_fcs = glob.glob(glob_fcs)

    if len(fn_fcs) == 0:
        raise IOError('FCS file not found')
    if len(fn_fcs) > 1:
        raise IOError('Multiple FCS files found')
    fn_fcs = fn_fcs[0]

    glob_index = fdn+'*'+plate+'*_Index.csv'
    fn_index = glob.glob(glob_index)
    if len(fn_index) == 0:
        raise IOError('Index file not found')
    if len(fn_index) > 1:
        raise IOError('Multiple index files found')
    fn_index = fn_index[0]

    meta, data = fcsparser.parse(fn_fcs, reformat_meta=True)
    out['fcs_meta'] = meta
    out['fcs_data'] = data

    data_index = pd.read_csv(fn_index, sep=',', index_col='Index')
    data_index.index = pd.Index(
            [plate+'_'+i for i in data_index.index],
            name='name')
    out['index_data'] = data_index

    # Post-processing:
    # - index data may use BSC instead of SSC
    if 'BSC-A' in data_index.columns:
        data_index.rename(
                columns={
                    'BSC-A': 'SSC-A',
                    'BSC-H': 'SSC-H',
                    'BSC-W': 'SSC-W'},
                inplace=True)

    # - index data may use FITC-Compensated, rename to FITC
    rename_comp = {}
    for x in data_index.columns:
        if x.endswith('-Compensated'):
            rename_comp[x] = x[:-len('-Compensated')]
    data_index.rename(columns=rename_comp, inplace=True)

    # - index data is compensated, but FCS is not by default
    if '$COMP' in out['fcs_meta']:
        com = list(map(float, out['fcs_meta']['$COMP'].split(',')))[1:]
        n_chan = int(np.sqrt(len(com)))
        # This is the right order, Fortran style I guess
        com = np.reshape(com, (n_chan, n_chan)).T
        fluors = out['fcs_meta']['_channel_names_'][-n_chan:]
        com = pd.DataFrame(
                data=com,
                index=fluors,
                columns=fluors)
        out['fcs_data'].loc[:, fluors] = out['fcs_data'].loc[:, fluors].values @ com.values

    # - add plate name to all FCS lines
    out['fcs_data'].loc[:, 'plate'] = plate

    return out


def parse_annotations(tissue):
    import glob
    if 'annotation glob' in config[tissue]:
        glb = config[tissue]['annotation glob']
    else:
        glb = tissue
    fn_glb = '../../data/MACAtSNE/{:}FACSmap.csv'.format(glb)

    fns = glob.glob(fn_glb)
    if len(fns) == 0:
        raise IOError('Annotation file not found for tissue: {:}'.format(tissue))
    elif len(fns) > 1:
        raise IOError('Several annotation files found for tissue: {:}'.format(tissue))
    else:
        fn = fns[0]

    out = pd.read_csv(fn, sep=',', index_col=0)

    #import ipdb; ipdb.set_trace()

    # Rename columns
    # Sometimes there is no subannotation
    out.rename(columns={
        'cluster.annotation': 'annotation',
        'cluster.subannotation': 'subannotation',
        'cluster': 'annotation'},
        inplace=True)

    if 'subannotation' not in out.columns:
        out['subannotation'] = np.nan

    ctc = []
    for key, datum in out.iterrows():
        ct = datum['annotation']
        cst = str(datum['subannotation'])
        if (cst != ct) and (cst.lower() != 'nan'):
            ct += ', '+cst
        ctc.append(ct)
    out.loc[:, 'cell_type_call'] = ctc

    # NOTE: some typos around
    out.loc[out['cell_type_call'] == 'Immunue cells', 'cell_type_call'] = 'Immune cells'

    # NOTE: brain plate 937 is probably a misannotation for brain plate 585
    out.index = [i if 'MAA000937' not in i else i.split('.')[0]+'.'+'MAA000585'
            for i in out.index]
    out['plate'] = [i.split('.')[1] for i in out.index]

    out.index.name = 'id'
    out['name'] = ['{1}_{0}'.format(*i.split('.')) for i in out.index]

    out.set_index('name', drop=True, inplace=True)

    plates = np.unique(out['plate'])

    return out, plates


def parse_counts(tissue):
    import glob
    if 'annotation glob' in config[tissue]:
        glb = config[tissue]['annotation glob']
    else:
        glb = tissue
    fn_glb = '../../data/MACAtSNE/{:}CountTable.csv'.format(glb)
    fns = glob.glob(fn_glb)
    if len(fns) == 0:
        raise IOError('Annotation file not found for tissue: {:}'.format(tissue))
    elif len(fns) > 1:
        raise IOError('Several annotation files found for tissue: {:}'.format(tissue))
    else:
        fn = fns[0]

    out = pd.read_csv(fn, sep=',', index_col=0)
    out.columns = ['{1}_{0}'.format(*(c.split('.')[:2])) for c in out.columns]
    out.index.name = 'GeneName'
    out.columns.name = 'Cell'
    return out


def parse_go_plasma_membrane():
    # GeneNames are unique, I checked
    fn = '../../data/go/plasma_membrane.tsv'
    out = pd.read_csv(fn, sep='\t', usecols=[0, 1], index_col=1).index
    return out


def get_dataset(tissue, membrane_only=True):
    counts = parse_counts(tissue)
    if membrane_only:
        go = parse_go_plasma_membrane()
        genes_membrane = go[go.isin(counts.index)]
        counts = counts.loc[genes_membrane]

    ds = Dataset(
            samplesheet=SampleSheet(cell_types),
            counts_table=CountsTable(counts),
            )
    return ds


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('tissues', nargs='+',
                        help='tissues to study')
    parser.add_argument('--save', action='store_true',
                        help='Store to file instead of showing')
    parser.add_argument('--explore', action='store_true',
                        help='Just explore the FACS data, do not analyze')
    parser.add_argument('--maxplates', type=int, default=0,
                        help='Max number of plates to plot')
    parser.add_argument('--groupby', nargs='+', default=('cell type call',),
                        choices=('cell type call', 'Good quality'),
                        help='Group cells by these criteria')
    parser.add_argument('--with-comparison', action='store_true',
                        help='Compare surface protein with transcript levels')
    args = parser.parse_args()

    if (len(args.tissues) == 1) and (args.tissues[0] == 'all'):
        args.tissues = tuple(config.keys())

    plate_meta = parse_plate_metadata()
    sys.exit()

    markers = {}
    for tissue in args.tissues:
        print(tissue)
        cell_types, plates = parse_annotations(tissue)

        if not args.explore:
            print('{:} has plates:'.format(tissue))
            print('\n'.join(plates))

        facs_data = {}
        for ip, plate in enumerate(plates):
            print('n. {:}, {:}'.format(ip+1, plate))

            # Some plates are missing
            try:
                facs_datum = parse_facs_plate(tissue, plate)
            except IOError:
                print('Missing')
                continue

            if args.explore:
                print('Found')
                continue

            ind_bool = facs_datum['index_data'].index.isin(cell_types.index)
            ind = facs_datum['index_data'].index[ind_bool]
            facs_datum['index_data']['Good quality'] = ind_bool

            facs_datum['index_data']["cell type call"] = "missing"
            facs_datum['index_data'].loc[ind, "cell type call"] = \
                cell_types.loc[ind, "cell_type_call"]

            if 'Reads' in cell_types.columns:
                facs_datum['index_data']["n reads"] = 0
                facs_datum['index_data'].loc[ind, "n reads"] = cell_types.loc[ind, "Reads"]

            if 'Genes' in cell_types.columns:
                facs_datum['index_data']["n genes"] = 0
                facs_datum['index_data'].loc[ind, "n genes"] = cell_types.loc[ind, "Genes"]

            if 'tSNEx' in cell_types.columns:
                facs_datum['index_data']["tSNEx"] = 0.0
                facs_datum['index_data'].loc[ind, "tSNEx"] = cell_types.loc[ind, "tSNEx"]
                facs_datum['index_data']["tSNEy"] = 0.0
                facs_datum['index_data'].loc[ind, "tSNEy"] = cell_types.loc[ind, "tSNEy"]

            facs_data[plate] = facs_datum

        if args.explore:
            continue

        cell_types_total = np.unique(cell_types['cell_type_call'])
        cell_types_total = np.append(cell_types_total, ['missing'])
        n_colors = len(cell_types_total)
        colors = {
                'Good quality': {
                    True: 'green',
                    False: 'red'},
                'cell type call': dict(zip(
                    cell_types_total,
                    sns.color_palette('husl', n_colors=n_colors))),
                }
        dead_stain = config[tissue]['dead stain']

        if args.with_comparison:
            ds = get_dataset(tissue, membrane_only=True)
            ds.counts.normalize(inplace=True)

        ip = 0
        for plate in plates:
            if plate not in facs_data:
                continue
            ip += 1
            facs_datum = facs_data[plate]
            for ig, groupby in enumerate(args.groupby):
                n_groups = len(np.unique(facs_datum['index_data'][groupby]))

                nplots = 2
                if ('dead stain' in config[tissue]) and (config[tissue]['dead stain'] is not None):
                    has_dead = True
                    nplots += 1
                else:
                    has_dead = False

                if 'plots' in config[tissue]:
                    nplots += len(config[tissue]['plots'])

                    if args.with_comparison:
                        ch_comp = []
                        ab_comp = []
                        gene_comp = []
                        for pl in config[tissue]['plots']:
                            for chname in (pl['x'], pl['y']):
                                # Annotated axes
                                if ':' in chname:
                                    # Multiple antibodies in the same channel are useless
                                    if '/' in chname:
                                        continue
                                    # Keep consistent case for antibodies
                                    abname = chname.split(':')[0].upper()
                                    gname = channels_to_genes.get(abname, abname)
                                    if gname not in ds.counts.index:
                                        continue
                                    ch_comp.append(chname)
                                    ab_comp.append(abname)
                                    gene_comp.append(gname)
                        print(ab_comp)
                        nplots += len(ab_comp)

                if nplots == 2:
                    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
                elif nplots == 3:
                    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(17, 4))
                elif nplots == 4:
                    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 9))
                    axs = axs.ravel()
                elif nplots <= 6:
                    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(17, 9))
                    axs = axs.ravel()
                else:
                    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(17, 12))
                    axs = axs.ravel()

                ax = axs[0]
                cells_plotted = []
                for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                    if qual == 'missing':
                        continue
                    cells_plotted.extend(list(datum.index))
                    datum.plot(
                            kind='scatter',
                            x='tSNEx',
                            y='tSNEy',
                            ax=ax,
                            color=colors[groupby][qual],
                            zorder=5,
                            alpha=0.7,
                            )
                other_cells = cell_types.loc[~cell_types.index.isin(cells_plotted)]
                other_cells.plot(
                        kind='scatter',
                        x='tSNEx',
                        y='tSNEy',
                        ax=ax,
                        marker='s',
                        color=[0.4] * 3,
                        grid=True,
                        zorder=4,
                        alpha=0.07,
                        s=10,
                        )
                ax.grid(True)

                ax = axs[1]
                posx = 'FSC-A'
                posy = 'SSC-A'
                if 'xlim' in config[tissue] and posx in config[tissue]['xlim']:
                    xlim = config[tissue]['xlim'][posx]
                else:
                    xlim = (1e2, 1e6)
                if 'ylim' in config[tissue] and posy in config[tissue]['ylim']:
                    ylim = config[tissue]['ylim'][posy]
                else:
                    ylim = (1e3, 1e6)

                for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                    if qual == 'missing':
                        continue
                    datum.plot(
                            kind='scatter',
                            x=posx,
                            y=posy,
                            label=qual,
                            ax=ax,
                            color=colors[groupby][qual],
                            zorder=5,
                            alpha=0.5,
                            )

                facs_datum['fcs_data'].plot(
                        kind='scatter',
                        x='FSC-A',
                        y='SSC-A',
                        ax=ax,
                        logy=True,
                        xlim=xlim,
                        ylim=ylim,
                        color=[0.4] * 3,
                        grid=True,
                        zorder=4,
                        alpha=0.01,
                        s=10,
                        )

                legend_kwargs = {
                        'loc': 'lower right',
                        'title': groupby,
                }
                if 'legend' in config[tissue]:
                    legend_kwargs.update(config[tissue]['legend'])
                ax.legend(**legend_kwargs)

                if has_dead:
                    ax = axs[2]
                    posy = dead_stain
                    if 'ylim' in config[tissue] and posy in config[tissue]['ylim']:
                        ylim = config[tissue]['ylim'][posy]
                    else:
                        ylim = (1e1, 1e6)

                    for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                        if qual == 'missing':
                            continue
                        datum.plot(
                                kind='scatter',
                                x='TIME',
                                y=dead_stain,
                                ax=ax,
                                color=colors[groupby][qual],
                                zorder=5,
                                alpha=0.7,
                                )

                    facs_datum['fcs_data'].plot(
                            kind='scatter',
                            x='TIME',
                            y=dead_stain,
                            ax=ax,
                            logy=True,
                            ylim=ylim,
                            color=[0.4] * 3,
                            grid=True,
                            zorder=4,
                            alpha=0.02,
                            s=10,
                            )

                if 'plots' in config[tissue]:
                    for ipl, pdict in enumerate(config[tissue]['plots']):
                        posx = pdict['x']
                        posy = pdict['y']
                        if 'xlim' in config[tissue] and posx in config[tissue]['xlim']:
                            xlim = config[tissue]['xlim'][posx]
                        else:
                            xlim = (1e2, 1e5)
                        if 'ylim' in config[tissue] and posy in config[tissue]['ylim']:
                            ylim = config[tissue]['ylim'][posy]
                        else:
                            ylim = (1e1, 1e6)

                        ax = axs[ipl + int(has_dead) + 2]
                        for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                            if qual == 'missing':
                                continue
                            datum.plot(
                                    kind='scatter',
                                    x=posx,
                                    y=posy,
                                    ax=ax,
                                    color=colors[groupby][qual],
                                    zorder=5,
                                    alpha=0.7,
                                    )

                        facs_datum['fcs_data'].plot(
                                kind='scatter',
                                x=posx,
                                y=posy,
                                ax=ax,
                                logx=(posx not in ('TIME', 'FSC-A')),
                                logy=(posy not in ('TIME', 'FSC-A')),
                                xlim=xlim,
                                ylim=ylim,
                                color='grey',
                                edgecolor='none',
                                grid=True,
                                zorder=4,
                                alpha=0.02,
                                s=10,
                                )

                        if 'antibodies' in config[tissue]:
                            if posx in config[tissue]['antibodies']:
                                ax.set_xlabel(config[tissue]['antibodies'][posx])
                            if posy in config[tissue]['antibodies']:
                                ax.set_ylabel(config[tissue]['antibodies'][posy])

                if args.with_comparison:
                    for ipl, (chname, abname, gname) in enumerate(zip(ch_comp, ab_comp, gene_comp)):
                        ind = np.intersect1d(
                                facs_datum['index_data'].index,
                                ds.counts.columns,
                                )
                        x = facs_datum['index_data'].loc[ind, chname]
                        y = ds.counts.loc[gname, ind] + ds.counts.pseudocount
                        rho = spearmanr(x, y)[0]
                        if ('xlim' in config[tissue]) and (chname in config[tissue]['xlim']):
                            xlim = config[tissue]['xlim'][chname]
                        else:
                            xlim = (1e1, 1e6)
                        ax = axs[ipl + int(has_dead) + 2 + len(config[tissue]['plots'])]
                        ax.scatter(
                                x, y,
                                s=20,
                                color='steelblue',
                                label='$\\rho = {:.2f}$'.format(rho),
                                )
                        ax.legend(loc='best', fontsize='8')
                        ax.set_xlabel(abname+', FACS stain')
                        ax.set_ylabel(gname+', transcript level')
                        ax.grid(True)
                        ax.set_xlim(*xlim)
                        ax.set_ylim(ds.counts.pseudocount * 0.9, 2e5)
                        ax.set_xscale('log')
                        ax.set_yscale('log')

                subtissue = plate_meta.loc[plate, 'subtissue']
                if subtissue:
                    title = '{:}, {:}, {:}'.format(plate, tissue, subtissue)
                else:
                    title = '{:}, {:}'.format(plate, tissue)
                fig.suptitle(title)
                plt.tight_layout(rect=(0, 0, 1, 0.96))

                if args.save:
                    fig.savefig(
                        '../../figures/second_screen/{:}_{:}_{:}.png'.format(tissue, plate, str(groupby)))
                    plt.close(fig)

            if ip == args.maxplates:
                break

    if not args.save:
        plt.ion()
        plt.show()
