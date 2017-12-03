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
    if 'brain' in tissue:
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
    parser.add_argument('--regenerate', action='store_true',
                        help='Store to file instead of showing')
    parser.add_argument('--save', action='store_true',
                        help='Store to file instead of showing')
    parser.add_argument('--maxplates', type=int, default=0,
                        help='Max number of plates to plot')
    args = parser.parse_args()

    if (len(args.tissues) == 1) and (args.tissues[0] == 'all'):
        args.tissues = tuple(config.keys())

    plate_meta = parse_plate_metadata()

    fn_cache = '../../data/correlations_ab_transcript.tsv'
    if not os.path.isfile(fn_cache) or args.regenerate:
        correlations = []
        for tissue in args.tissues:
            print(tissue)
            print('Load annotations')
            cell_types, plates = parse_annotations(tissue)

            print('Load counts')
            ds = get_dataset(tissue, membrane_only=True)

            print('Checking stains for comparison')
            if 'plots' not in config[tissue]:
                continue
            ch_comp = []
            ab_comp = []
            gene_comp = []
            for pl in config[tissue]['plots']:
                for chname in (pl['x'], pl['y']):
                    if ('antibodies' in config[tissue]) and (chname in config[tissue]['antibodies']):
                        chfull = config[tissue]['antibodies'][chname]
                    else:
                        chfull = chname

                    print('Checking {:}'.format(chfull))

                    # Annotated axes
                    if ':' in chfull:
                        # Multiple antibodies in the same channel are useless
                        if '/' in chfull:
                            continue
                        # Keep consistent case for antibodies
                        abname = chfull.split(':')[0]
                        if abname.startswith('Cd'):
                            abname = 'CD'+abname[2:]
                        gname = channels_to_genes.get(abname, abname)
                        if gname not in ds.counts.index:
                            continue
                        if gname in gene_comp:
                            continue
                        ch_comp.append(chname)
                        ab_comp.append(abname)
                        gene_comp.append(gname)

            if not len(gene_comp):
                continue

            print('Normalize counts')
            ds.counts.normalize(inplace=True)


            facs_data = {}
            ip = 0
            for plate in plates:
                ip += 1
                subtissue = plate_meta.loc[plate, 'subtissue']

                print('{:}, {:}'.format(plate, subtissue))

                # Some plates are missing
                try:
                    facs_datum = parse_facs_plate(tissue, plate)
                except IOError:
                    print('Missing')
                    continue

                ind_bool = facs_datum['index_data'].index.isin(cell_types.index)
                ind = facs_datum['index_data'].index[ind_bool]
                facs_datum['index_data']['Good quality'] = ind_bool

                facs_datum['index_data']["cell type call"] = "missing"
                facs_datum['index_data'].loc[ind, "cell type call"] = \
                    cell_types.loc[ind, "cell_type_call"]

                facs_data[plate] = facs_datum

                facs_datum = facs_data[plate]

                for ipl, (chname, abname, gname) in enumerate(zip(ch_comp, ab_comp, gene_comp)):
                    ind = np.intersect1d(
                            facs_datum['index_data'].index,
                            ds.counts.columns,
                            )
                    x = facs_datum['index_data'].loc[ind, chname]
                    y = ds.counts.loc[gname, ind] + ds.counts.pseudocount
                    ctype = facs_datum['index_data'].loc[ind, 'cell type call']
                    df = pd.concat([x, y, ctype], axis=1)
                    rho = spearmanr(x, y)[0]
                    xl = np.log10(0.1 + x.values)
                    yl = np.log10(0.1 + y.values)
                    n_cells = len(xl)
                    correlations.append({
                        'rho': rho,
                        'n_cells': n_cells,
                        'tissue': tissue,
                        'subtissue': subtissue,
                        'plate': plate,
                        'gene': gname,
                        'ab_mean': xl.mean(),
                        'ab_std': xl.std(),
                        'gene_mean': yl.mean(),
                        'gene_std': yl.std(),
                        })

                if ip == args.maxplates:
                    break

        # NOTE: both axes got log10-ed before calculating mean and std
        correlations = pd.DataFrame(correlations)

    else:
        correlations = pd.read_csv(fn_cache, sep='\t', index_col=0)

    # FIXME: why are there NaNs at all?
    correlations.dropna(subset=('ab_mean', 'ab_std', 'rho'), inplace=True)

    fig, ax = plt.subplots(figsize=(7, 5))

    genes = np.unique(correlations['gene'])
    colors = sns.color_palette('husl', n_colors=len(genes))
    for ig, (gene, corr) in enumerate(correlations.groupby('gene')):
        corr.plot(
                kind='scatter',
                x='ab_std',
                y='rho',
                ax=ax,
                color=colors[ig],
                label=gene,
                )
    ax.grid(True)
    ax.set_ylim(-1, 1)
    ax.set_ylabel('Spearman rho')
    ax.set_xlabel('Antibody stain\nstandard dev (decades)')

    x = correlations['ab_std'].values
    y = correlations['rho'].values
    m = (x @ y) / (x @ x)
    xfit = np.linspace(0, 1, 100)
    ax.plot(xfit, m * xfit, lw=2, color='k', label='Fit, m={:.2f}'.format(m))
    ax.legend(loc='best')

    plt.tight_layout(rect=(0, 0.02, 1, 1))

    if args.save:
        fig.savefig('../../figures/correlations_vs_abstd.png')
        plt.close(fig)
    else:
        plt.ion()
        plt.show()
