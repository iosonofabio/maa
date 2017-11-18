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
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset
from singlet.counts_table import CountsTable
from singlet.samplesheet import SampleSheet


# Globals
config = {
    'colon': {
        'dead stain': 'CD45-DAPI: Pacific Blue-A',
        'plots': [{
            'x': 'CD44: APC-A',
            'y': 'CD66a: PE-A',
            }, {
            'x': 'CD44: APC-A',
            'y': '',
            }],
        'xlim': {'CD44: APC-A': (2e2, 5e4)},
        'ylim': {'CD66a: PE-A': (7e2, 5e5)},
        },
    'trachea': {
        'dead stain': 'dead: Sytox blue-A',
        },
    'pancreas': {
        'dead stain': 'dead: APC-Cy7-A',
        'legend': {
            'fontsize': 8,
            },
        },
    'heart': {
        'dead stain': 'PI-A',
        'antibodies': {
            'PI-A': 'dead: PI-A',
            },
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        'xlim': {'FSC-A': (1, 8e5)},
        'ylim': {'SSC-A': (5e3, 1e6)},
        },
    'aorta': {
        'dead stain': 'PI-A',
        'antibodies': {
            'PI-A': 'dead: PI-A',
            },
        },
    'spleen': {
        'annotation glob': 'Spleen',
        'dead stain': 'PI-A',
        'antibodies': {
            'PI-A': 'dead: PI-A',
            },
        'xlim': {'FSC-A': (1, 9.5e5)},
        'ylim': {'SSC-A': (5e3, 1e6)},
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        },
    'tongue': {
        'dead stain': 'Brilliant Violet 421-A',
        'plots': [{
            'x': 'FITC-A',
            'y': 'APC-A',
            }],
        'antibodies': {
            'Brilliant Violet 421-A': 'dead: Brilliant Violet 421-A',
            # FIXME: ask them this thing!
            'FITC-A': '???: FITC-A',
            'APC-A': '???: APC-A',
            },
        'xlim': {'FITC-A': (3e2, 1e5)},
        'ylim': {'APC-A': (1e2, 1e5)},
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        },
    'bladder': {
        'annotation glob': 'Bladder',
        'dead stain': 'Brilliant Violet 421-A',
        'plots': [{
            'x': 'FITC-A',
            'y': 'APC-A',
            }],
        'antibodies': {
            'Brilliant Violet 421-A': 'dead: Brilliant Violet 421-A',
            # FIXME: ask them this thing!
            'FITC-A': '???: FITC-A',
            'APC-A': '???: APC-A',
            },
        'legend': {
            'ncol': 2,
            },
        },
    'brain_neuron': {
        'annotation glob': 'BrainNeuron',
        'dead stain': 'Lineage: Brilliant Violet 421-A',
        },
    'brain_microglia': {
        'annotation glob': 'brainMicroglia',
        'dead stain': 'Live/dead: PI-A',
        'plots': [{
            'x': 'CD11b: Brilliant Violet 421-A',
            'y': 'CD45: PE-Cy7-A',
            }],
        },
    'kidney': {
        'dead stain': 'dead: PI-A',
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        },
    'skin': {
        'dead stain': 'Pacific Blue-A',
        'plots': [{
            'x': 'FITC-A',
            'y': 'APC-A',
            }],
        'antibodies': {
            'Pacific Blue-A': 'dead: Pacific Blue-A',
            'FITC-A': 'a6 Integrin: FITC-A',
            'APC-A': 'CD34: Alexa647-A',
            },
        'ylim': {'APC-A': (1e1, 7e3)},
        },
    'fat': {
        'dead stain': 'Pacific Blue-A',
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        'plots': [{
            'x': 'FITC-A',
            'y': 'PE-Cy7-A',
            }, {
            'x': 'APC-A',
            'y': 'FSC-A',
            }],
        'antibodies': {
            'Pacific Blue-A': 'dead: Pacific Blue-A',
            'FITC-A': 'CD31: FITC-A',
            'APC-A': 'SCA-1: APC-A',
            'PE-Cy7-A': 'CD45: PE-Cy7-A',
            },
        },
    'liver': {
        },
    'lung': {
        },
    # Muscle is all ARIA sorting at the VA
    'diaphragm': {
        'sorter': 'ARIA',
        'dead stain': None,  # ??
        'antibodies': {
            'Pacific Blue-A': 'SCA-1: Pacific Blue-A',
            'FITC-A': 'CD31: FITC-A',
            'APC-A': 'CD45: APC-A',
            'PE-Cy7-A': 'VCAM: PE-Cy7-A',
            },
        'xlim': {'FSC-A': (1e2, 2.7e5)},
        'ylim': {'SSC-A': (3e3, 3e5)},
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        },
    'muscle': {
        'sorter': 'ARIA',
        'dead stain': None,  # ??
        'antibodies': {
            'Pacific Blue-A': 'Ly-6A/E: Pacific Blue-A',
            'FITC-A': 'CD45: FITC-A',
            'APC-A': 'CD31: APC-A',
            'PE-Cy7-A': 'CD106: PE-Cy7-A',
            },
        'xlim': {'FSC-A': (1e2, 2.67e5)},
        'ylim': {'SSC-A': (3e3, 3e5)},
        'legend': {
            'ncol': 2,
            'fontsize': 8,
            },
        },
    }


# Functions
def parse_facs_plate(tissue, plate):
    sorter = config[tissue].get('sorter', 'Sony')
    if sorter == 'Sony':
        return parse_facs_plate_sony(plate)
    else:
        return parse_facs_plate_aria(plate)


def parse_facs_plate_aria(plate):
    import glob
    import fcsparser

    out = {}

    fdn = '../../data/MACAFACS/index_fcs_*/'
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


def parse_facs_plate_sony(plate):
    import glob
    import fcsparser

    out = {}

    fdn = '../../data/MACAFACS/index_fcs_*/'
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
    parser.add_argument('--find-markers', action='store_true',
                        help='Find new markers, do not plot')
    parser.add_argument('--maxplates', type=int, default=0,
                        help='Max number of plates to plot')
    parser.add_argument('--groupby', nargs='+', default=('cell type call',),
                        choices=('cell type call', 'Good quality'),
                        help='Group cells by these criteria')
    args = parser.parse_args()

    markers = {}
    for tissue in args.tissues:
        print(tissue)
        cell_types, plates = parse_annotations(tissue)

        if not args.explore:
            ds = get_dataset(tissue, membrane_only=True)
            ds.counts.normalize(inplace=True)

        if args.find_markers:
            ds.counts.log(base=10, inplace=True)
            for ct in np.unique(cell_types['cell_type_call']):
                cn = 'cell_type_{:}'.format(ct)
                ds.samplesheet[cn] = False
                ds.samplesheet.loc[ds.samplesheet['cell_type_call'] == ct, cn] = True

            for ct in np.unique(cell_types['cell_type_call']):
                cn = 'cell_type_{:}'.format(ct)
                dst = ds.split(cn)
                res = dst[True].compare(dst[False])
                med = dst[True].counts.median(axis=1)
                res['median diff in log10'] = med - dst[False].counts.median(axis=1)
                res['median exp in subtype'] = med
                qua = dst[True].counts.quantile([0.1, 0.9], axis=1).T
                res['10% exp in subtype'] = qua[0.1]
                res['90% exp in subtype'] = qua[0.9]
                ma_ind = (res['P-value'] < 1e-6) & (res['median exp in subtype'] >= 2)
                ma = res.loc[ma_ind].sort_values(by='P-value').index[:10]
                markers[(tissue, ct, '+')] = ', '.join(ma)

            continue

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

                if nplots == 2:
                    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
                elif nplots == 3:
                    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(17, 4))
                elif nplots == 4:
                    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 9))
                    axs = axs.ravel()
                else:
                    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(17, 9))
                    axs = axs.ravel()

                ax = axs[0]
                for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                    if qual == 'missing':
                        continue
                    datum.plot(
                            kind='scatter',
                            x='tSNEx',
                            y='tSNEy',
                            ax=ax,
                            color=colors[groupby][qual],
                            zorder=5,
                            alpha=0.7,
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
                                logx=True,
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

                fig.suptitle(plate+', {:}'.format(tissue))
                plt.tight_layout(rect=(0, 0, 1, 0.96))

                if args.save:
                    fig.savefig(
                        '../../figures/first_screen/{:}_{:}_{:}.png'.format(tissue, plate, str(groupby)))
                    plt.close(fig)

            if ip == args.maxplates:
                break

    if args.find_markers:
        markers = pd.Series(markers)
        markers.index.names = ['tissue', 'cell type', 'positive/negative stain']

    if not args.save:
        plt.ion()
        plt.show()
