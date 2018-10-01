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
    md = md.reset_index().dropna(subset=['name']).set_index('name')

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


def parse_facs_plate(plate, maxfcs=0, dry=False):
    import glob

    # Check Aria
    glob_index = '../../data/MACAFACS/index_fcs_*/'+plate+'_Index.fcs'
    fn_index = glob.glob(glob_index)
    if len(fn_index) > 1:
        raise IOError('Multiple index files found')
    elif len(fn_index) == 0:
        # Check Sony
        glob_index = '../../data/MACAFACS/index_fcs_*/*'+plate+'*_Index.csv'
        fn_index = glob.glob(glob_index)
        if len(fn_index) > 1:
            raise IOError('Multiple index files found')
        elif len(fn_index) == 0:
            raise IOError('Index file not found')
        else:
            sorter = 'Sony'
    else:
        sorter = 'Aria'

    if dry:
        print('Plate {:}, sorter {:}'.format(plate, sorter))
        return

    if sorter == 'Sony':
        return parse_facs_plate_sony(plate, glb='*', maxfcs=maxfcs)
    else:
        return parse_facs_plate_aria(plate, maxfcs=maxfcs)


def parse_facs_plate_aria(plate, glb='*', maxfcs=50000):
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

    if maxfcs != 0:
        meta, data = fcsparser.parse(fn_fcs, reformat_meta=True)
        if maxfcs != -1:
            data = data.loc[:maxfcs]
        out['fcs_meta'] = meta
        out['fcs_data'] = data

    glob_index = fdn+plate+'_Index.fcs'
    fn_index = glob.glob(glob_index)
    if len(fn_index) == 0:
        raise IOError('Index file not found')
    if len(fn_index) > 1:
        raise IOError('Multiple index files found')
    fn_index = fn_index[0]

    # Index data
    meta_index, data_index = fcsparser.parse(
            fn_index,
            meta_data_only=False,
            reformat_meta=True)
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


def parse_facs_plate_sony(plate, glb='*', maxfcs=50000):
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

    if maxfcs != 0:
        meta, data = fcsparser.parse(fn_fcs, reformat_meta=True)
        if maxfcs != -1:
            data = data.loc[:maxfcs]
        # add plate name to all FCS lines
        data.loc[:, 'plate'] = plate
        # compensate FCS data
        if '$COMP' in meta:
            com = list(map(float, meta['$COMP'].split(',')))[1:]
            n_chan = int(np.sqrt(len(com)))
            # This is the right order, Fortran style I guess
            com = np.reshape(com, (n_chan, n_chan)).T
            fluors = meta['_channel_names_'][-n_chan:]
            com = pd.DataFrame(
                    data=com,
                    index=fluors,
                    columns=fluors)
            data.loc[:, fluors] = data.loc[:, fluors].values @ com.values
        out['fcs_meta'] = meta
        out['fcs_data'] = data

    # Index data
    glob_index = fdn+'*'+plate+'*_Index.csv'
    fn_index = glob.glob(glob_index)
    if len(fn_index) == 0:
        raise IOError('Index file not found')
    if len(fn_index) > 1:
        raise IOError('Multiple index files found')
    fn_index = fn_index[0]

    data_index = pd.read_csv(fn_index, sep=',', index_col='Index')
    data_index.index = pd.Index(
            [plate+'_'+i for i in data_index.index],
            name='name')
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
    data_index.loc[:, 'plate'] = plate
    out['index_data'] = data_index

    return out


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    g = parser.add_mutually_exclusive_group()
    g.add_argument('--dry', action='store_true',
                   help='Dry run (just check file existance)')
    g.add_argument('--save', action='store_true',
                   help='Store to file instead of showing')
    args = parser.parse_args()

    plate_meta = parse_plate_metadata()

    for plate, metadatum in plate_meta.iterrows():
        tissue = metadatum['tissue']
        try:
            datum = parse_facs_plate(plate, maxfcs=0, dry=args.dry)
        except IOError:
            continue

        if args.save:
            datum = datum['index_data']
            for col in ['date.sorted', 'tissue', 'mouse.id',
                        'FACS.selection', 'nozzle.size',
                        'FACS.instument', 'mouse.age',
                        'n.wells']:
                datum[col] = metadatum[col]

            fn_out = '../../data/MACAFACS_index_tsv/{:}.tsv'.format(plate)
            print('Plate {:}: save to CSV'.format(plate))
            datum.to_csv(
                    fn_out,
                    sep='\t',
                    index=True,
                    )

    #for tissue, subtissue in tissues:
    #    cell_types, plates = parse_annotations(tissue, subtissue=subtissue)

    #    if len(plates) == 0:
    #        raise ValueError('No plates found for {:}, {:}'.format(tissue, subtissue))

    #    print('{:}, {:} has plates:'.format(tissue, subtissue))
    #    print('\n'.join(plates))

    #    facs_data = {}
    #    for ip, plate in enumerate(plates):
    #        print('n. {:}, {:}'.format(ip+1, plate))

    #        # Some plates are missing
    #        try:
    #            facs_datum = parse_facs_plate(tissue, plate, maxfcs=args.maxfcs)
    #        except IOError:
    #            print('Missing')
    #            continue

    #        ind_bool = facs_datum['index_data'].index.isin(cell_types.index)
    #        ind = facs_datum['index_data'].index[ind_bool]
    #        facs_datum['index_data']['Good quality'] = ind_bool

    #        facs_datum['index_data']["cell type call"] = "missing"
    #        facs_datum['index_data'].loc[ind, "cell type call"] = \
    #            cell_types.loc[ind, "cell_type_call"]

    #        facs_datum['index_data']['all_true'] = True

    #        if 'Reads' in cell_types.columns:
    #            facs_datum['index_data']["n reads"] = 0
    #            facs_datum['index_data'].loc[ind, "n reads"] = cell_types.loc[ind, "Reads"]

    #        if 'Genes' in cell_types.columns:
    #            facs_datum['index_data']["n genes"] = 0
    #            facs_datum['index_data'].loc[ind, "n genes"] = cell_types.loc[ind, "Genes"]

    #        if 'tSNEx' in cell_types.columns:
    #            facs_datum['index_data']["tSNEx"] = 0.0
    #            facs_datum['index_data'].loc[ind, "tSNEx"] = cell_types.loc[ind, "tSNEx"]
    #            facs_datum['index_data']["tSNEy"] = 0.0
    #            facs_datum['index_data'].loc[ind, "tSNEy"] = cell_types.loc[ind, "tSNEy"]

    #        facs_data[plate] = facs_datum

    #    if args.mergeplates:
    #        tmp_index = pd.concat([facs_datum['index_data'] for _, facs_datum in facs_data.items()], axis=0)
    #        if args.maxfcs != 0:
    #            tmp_fcs = pd.concat([facs_datum['fcs_data'] for _, facs_datum in facs_data.items()], axis=0)
    #        facs_data = {'all plates': {
    #            'index_data': tmp_index},
    #            }
    #        if args.maxfcs != 0:
    #            facs_data['fcs_data'] = tmp_fcs

    #    cell_types_total = np.unique(cell_types['cell_type_call'])
    #    cell_types_total = np.append(cell_types_total, ['missing'])
    #    n_colors = len(cell_types_total)
    #    colors = {
    #            'Good quality': {
    #                True: 'green',
    #                False: 'red'},
    #            'cell type call': dict(zip(
    #                cell_types_total,
    #                sns.color_palette('husl', n_colors=n_colors))),
    #            'all_true': {
    #                True: 'steelblue',
    #                }
    #            }
    #    dead_stain = config[tissue]['dead stain']

    #    ip = 0
    #    for plate in facs_data:
    #        ip += 1
    #        facs_datum = facs_data[plate]

    #        for ig, groupby in enumerate(args.groupby):
    #            n_groups = len(np.unique(facs_datum['index_data'][groupby]))

    #            nplots = 1
    #            if ('dead stain' in config[tissue]) and (config[tissue]['dead stain'] is not None):
    #                has_dead = True
    #                nplots += 1
    #            else:
    #                has_dead = False

    #            if 'plots' in config[tissue]:
    #                nplots += len(config[tissue]['plots'])

    #            if nplots == 1:
    #                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
    #                axs = [ax]
    #            elif nplots == 2:
    #                fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
    #            elif nplots == 3:
    #                fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(17, 4))
    #            elif nplots == 4:
    #                fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 9))
    #                axs = axs.ravel()
    #            elif nplots <= 6:
    #                fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(17, 9))
    #                axs = axs.ravel()
    #            elif nplots <= 9:
    #                fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(17, 12))
    #                axs = axs.ravel()
    #            else:
    #                fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(19, 12))
    #                axs = axs.ravel()

    #            ax = axs[0]
    #            posx = 'FSC-A'
    #            posy = 'SSC-A'
    #            if 'xlim' in config[tissue] and posx in config[tissue]['xlim']:
    #                xlim = config[tissue]['xlim'][posx]
    #            else:
    #                xlim = (1e2, 1e6)
    #            if 'ylim' in config[tissue] and posy in config[tissue]['ylim']:
    #                ylim = config[tissue]['ylim'][posy]
    #            else:
    #                ylim = (1e3, 1e6)

    #            for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
    #                if qual == 'missing':
    #                    continue
    #                datum.plot(
    #                        kind='scatter',
    #                        x=posx,
    #                        y=posy,
    #                        label=qual,
    #                        ax=ax,
    #                        color=colors[groupby][qual],
    #                        zorder=5,
    #                        alpha=0.5,
    #                        )

    #            if 'fcs_data' in facs_datum:
    #                facs_datum['fcs_data'].plot(
    #                        kind='scatter',
    #                        x='FSC-A',
    #                        y='SSC-A',
    #                        ax=ax,
    #                        color=[0.4] * 3,
    #                        zorder=4,
    #                        alpha=0.01,
    #                        s=10,
    #                        )
    #            ax.set_xlim(xlim)
    #            ax.set_ylim(ylim)
    #            ax.set_yscale('log')
    #            ax.grid(True)

    #            if args.groupby != ('all_true', ):
    #                legend_kwargs = {
    #                        'loc': 'lower right',
    #                        'title': groupby,
    #                }
    #                if 'legend' in config[tissue]:
    #                    legend_kwargs.update(config[tissue]['legend'])
    #                ax.legend(**legend_kwargs)
    #            else:
    #                ax.legend_.remove()

    #            if has_dead:
    #                ax = axs[1]
    #                posy = dead_stain
    #                if 'ylim' in config[tissue] and posy in config[tissue]['ylim']:
    #                    ylim = config[tissue]['ylim'][posy]
    #                else:
    #                    ylim = (1e1, 1e6)

    #                for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
    #                    if qual == 'missing':
    #                        continue
    #                    datum.plot(
    #                            kind='scatter',
    #                            x='TIME',
    #                            y=dead_stain,
    #                            ax=ax,
    #                            color=colors[groupby][qual],
    #                            zorder=5,
    #                            alpha=0.7,
    #                            )

    #                if 'fcs_data' in facs_datum:
    #                    facs_datum['fcs_data'].plot(
    #                            kind='scatter',
    #                            x='TIME',
    #                            y=dead_stain,
    #                            ax=ax,
    #                            color=[0.4] * 3,
    #                            zorder=4,
    #                            alpha=0.02,
    #                            s=10,
    #                            )
    #                ax.set_ylim(ylim)
    #                ax.set_yscale('log')
    #                ax.grid(True)

    #            if 'plots' in config[tissue]:
    #                for ipl, pdict in enumerate(config[tissue]['plots']):
    #                    posx = pdict['x']
    #                    posy = pdict['y']
    #                    if 'xlim' in config[tissue] and posx in config[tissue]['xlim']:
    #                        xlim = config[tissue]['xlim'][posx]
    #                    else:
    #                        xlim = (1e2, 1e5)
    #                    if 'ylim' in config[tissue] and posy in config[tissue]['ylim']:
    #                        ylim = config[tissue]['ylim'][posy]
    #                    else:
    #                        ylim = (1e1, 1e6)

    #                    ax = axs[ipl + int(has_dead) + 1]
    #                    for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
    #                        if qual == 'missing':
    #                            continue
    #                        datum.plot(
    #                                kind='scatter',
    #                                x=posx,
    #                                y=posy,
    #                                ax=ax,
    #                                color=colors[groupby][qual],
    #                                zorder=5,
    #                                alpha=0.7,
    #                                )

    #                    if 'fcs_data' in facs_datum:
    #                        facs_datum['fcs_data'].plot(
    #                                kind='scatter',
    #                                ax=ax,
    #                                x=posx,
    #                                y=posy,
    #                                color='grey',
    #                                edgecolor='none',
    #                                zorder=4,
    #                                alpha=0.02,
    #                                s=10,
    #                                )
    #                    ax.set_xlim(xlim)
    #                    ax.set_ylim(ylim)
    #                    if posx not in ('TIME', 'FSC-A'):
    #                        ax.set_xscale('log')
    #                    if posy not in ('TIME', 'FSC-A'):
    #                        ax.set_yscale('log')
    #                    ax.grid(True)

    #                    if 'antibodies' in config[tissue]:
    #                        if posx in config[tissue]['antibodies']:
    #                            ax.set_xlabel(config[tissue]['antibodies'][posx])
    #                        if posy in config[tissue]['antibodies']:
    #                            ax.set_ylabel(config[tissue]['antibodies'][posy])

    #            if subtissue is not None:
    #                title = '{:}, {:}, {:}'.format(plate, tissue, subtissue)
    #            else:
    #                title = '{:}, {:}'.format(plate, tissue)
    #            fig.suptitle(title)
    #            plt.tight_layout(rect=(0, 0, 1, 0.96))

    #            if args.save:
    #                fig.savefig(
    #                    '../../figures/facs_plots_maca/{:}_{:}_{:}.png'.format(tissue, plate, str(groupby)))
    #                plt.close(fig)

    #        if ip == args.maxplates:
    #            break

    #if not args.save:
    #    plt.ion()
    #    plt.show()
