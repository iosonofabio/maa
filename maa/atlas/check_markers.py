# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/17
content:    Check a priori markers for cell type separation.
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
    tissues_prediction = config['tissues_prediction']
    config = config['tissues']


# Functions
def parse_biolegend():
    fn = '../../data/ab_vendors/Biolegend.tsv'
    df = pd.read_csv(fn, sep='\t')

    if 'GeneName' not in df.columns:
        df.to_csv(fn+'.bak', sep='\t', index=False)

        from collections import Counter
        fn_conv = '../../data/ab_vendors/Biolegend_markers_conversion.tsv'
        df_conv = pd.read_csv(fn_conv, sep='\t')
        n_entries = Counter(df_conv['Input'].values)
        multiples = [k for (k, v) in n_entries.items() if v > 1]
        if len(multiples):
            print('Multiple entries:', multiples)
            raise ValueError('Some antibody target names have multiple entries')
        df_conv.set_index('Input', inplace=True)

        df['GeneName'] = ''
        newcols = df.columns[:2].tolist() + ['GeneName'] + df.columns[2:].tolist()
        df = df.loc[:, newcols]
        for k, datum in df.iterrows():
            if datum['Specificity'] in df_conv.index:
                df.loc[k, 'GeneName'] = df_conv.loc[datum['Specificity'], 'Symbol']

        df.to_csv(fn, sep='\t', index=False)
        print('New file saved to file')

    df.iloc[:, 3:] = (df.iloc[:, 3:] == '•')

    return df

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


def parse_counts(tissue, regenerate=False):
    import glob
    if 'annotation glob' in config[tissue]:
        glb = config[tissue]['annotation glob']
    else:
        glb = tissue

    if regenerate:
        cglbs = ('CountTable',)
    else:
        cglbs = ('CountTableNormalized', 'CountTable')

    for cglb in cglbs:
        fn_glb = '../../data/MACAtSNE/{:}{:}.csv'.format(glb, cglb)
        fns = glob.glob(fn_glb)
        if len(fns):
            break

    if len(fns) == 0:
        raise IOError('Counts file not found for tissue: {:}'.format(tissue))
    elif len(fns) > 1:
        raise IOError('Several counts files found for tissue: {:}'.format(tissue))
    else:
        fn = fns[0]

    out = pd.read_csv(fn, sep=',', index_col=0)
    if '.' in out.columns[0]:
        out.columns = ['{1}_{0}'.format(*(c.split('.')[:2])) for c in out.columns]
    out.index.name = 'GeneName'
    out.columns.name = 'Cell'
    out = CountsTable(out)
    if 'Normalized' in fn:
        out._normalized = 'counts_per_million'
    else:
        print('Normalize counts')
        out.normalize(inplace=True)
        print('Log counts')
        out.log(inplace=True)
        print('Write normalized counts to file')
        out.to_csv(fn[:-4]+'Normalized.csv', sep=',')
    return out


def parse_go_plasma_membrane():
    # GeneNames are unique, I checked
    fn = '../../data/go/plasma_membrane.tsv'
    out = pd.read_csv(fn, sep='\t', usecols=[0, 1, 3], index_col=1).iloc[:, :2]
    out.columns = ['GeneId', 'GONames']
    return out


def get_dataset(
        tissue,
        membrane_only=True,
        regenerate=False,
        go_contains=None,
        go_exclude=None):

    # Some tissues like brain were split for sorting, we merge them here
    dss = []
    for tissue_facs in tissues_prediction[tissue]:
        cell_types, plates = parse_annotations(tissue_facs)
        counts = parse_counts(tissue_facs, regenerate=regenerate)
        if membrane_only:
            go = parse_go_plasma_membrane().index
            genes_membrane = go[go.isin(counts.index)]
            counts = counts.loc[genes_membrane]

        if (go_contains is not None) and (go_exclude is not None):
            raise ValueError('Use either go_contains or go_exclude')
        if go_contains is not None:
            go = parse_go_plasma_membrane()
            genes = go.index[go['GONames'].str.contains(go_contains)]
            genes = np.intersect1d(genes, counts.index)
            counts = counts.loc[genes]
        elif go_exclude is not None:
            go = parse_go_plasma_membrane()
            genes = go.index[~go['GONames'].str.contains(go_exclude)]
            genes = np.intersect1d(genes, counts.index)
            counts = counts.loc[genes]

        dss.append({
            'samplesheet': cell_types,
            'counts': counts})

    if len(dss) == 1:
        ds = Dataset(
                samplesheet=SampleSheet(cell_types),
                counts_table=counts,
                )
        return ds
    else:
        # Merging is kind of messy because some genes are absent from either
        # subtissue (grrr); I put zeroes for now, Michelle is working on the
        # better solution (we have those numbers somewhere)
        genes = set()
        for ds in dss:
            genes |= set(ds['counts'].index.values)
        genes = pd.Index(sorted(genes), name=ds['counts'].index.name)
        for ds in dss:
            genes_missing = genes[~genes.isin(ds['counts'].index)]
            for gene in genes_missing:
                # The stuff is normalized, pseudocounted, and logged
                ds['counts'].loc[gene] = -1.0
            ds['counts'] = ds['counts'].loc[genes]
        ngenes = len(genes)
        ncells = sum(ds['samplesheet'].shape[0] for ds in dss)
        samplesheet_all = pd.concat([ds['samplesheet'] for ds in dss], axis=0)
        counts_all = pd.DataFrame(
                np.zeros((ngenes, ncells), float),
                index=genes,
                columns=samplesheet_all.index)
        for ds in dss:
            counts_all.loc[:, ds['counts'].columns.values] = ds['counts'].values
        counts_all = CountsTable(counts_all)
        if ds['counts']._normalized:
            counts_all._normalized = ds['counts']._normalized

        ds = Dataset(
                samplesheet=SampleSheet(samplesheet_all),
                counts_table=counts_all,
                )
        return ds


def make_subplots(nplots):
    if nplots == 1:
        fig, axs = plt.subplots(1, 1, figsize=(4, 4))
        axs = [axs]
    elif nplots == 2:
        fig, axs = plt.subplots(1, 2, figsize=(7, 4))
    elif nplots == 3:
        fig, axs = plt.subplots(1, 3, figsize=(9, 4))
    elif nplots == 4:
        fig, axs = plt.subplots(1, 4, figsize=(12, 4))
    elif nplots <= 6:
        fig, axs = plt.subplots(2, 3, figsize=(9, 7))
        axs = axs.ravel()
    elif nplots <= 8:
        fig, axs = plt.subplots(2, 4, figsize=(12, 7))
        axs = axs.ravel()
    elif nplots == 9:
        fig, axs = plt.subplots(3, 3, figsize=(9, 9))
        axs = axs.ravel()
    elif nplots <= 12:
        fig, axs = plt.subplots(3, 4, figsize=(12, 9))
        axs = axs.ravel()
    elif nplots <= 16:
        fig, axs = plt.subplots(4, 4, figsize=(12, 11))
        axs = axs.ravel()
    elif nplots <= 20:
        fig, axs = plt.subplots(4, 5, figsize=(12, 11))
        axs = axs.ravel()
    elif nplots <= 25:
        fig, axs = plt.subplots(5, 5, figsize=(12, 12))
        axs = axs.ravel()
    elif nplots <= 30:
        fig, axs = plt.subplots(5, 6, figsize=(14, 12))
        axs = axs.ravel()
    elif nplots <= 35:
        fig, axs = plt.subplots(5, 7, figsize=(15, 12))
        axs = axs.ravel()
    else:
        raise ValueError('Too many plots!')
    return fig, axs



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('tissues', nargs='+',
                        help='tissues to study')
    parser.add_argument('--regenerate', action='store_true',
                        help='Regenerate counts cache')
    parser.add_argument('--cell-types', nargs='+', default=None,
                        help='Limit to some cell types')
    parser.add_argument('--subtissue', default=None,
                        help='Limit to a subtissue. To split by subtissue use "all"')
    parser.add_argument('--markers', nargs='+', type=str,
                        help='Select markers')
    args = parser.parse_args()

    if (len(args.tissues) == 1) and (args.tissues[0] == 'all'):
        args.tissues = tuple(tissues_prediction.keys())

    if len(args.markers) > 2:
        raise ValueError('Only one or two markers supported')

    go = parse_go_plasma_membrane()

    # Get the list of commercially available antibodies
    ab_comm = []
    # Biolegend
    ab_comm_table = parse_biolegend()
    ab_unique = np.unique(ab_comm_table.dropna(subset=['GeneName'], axis=0)['GeneName'])
    ab_comm.append(ab_unique)
    # TODO: other vendors
    if len(ab_comm):
        ab_comm = np.unique(np.concatenate(ab_comm))

    plate_meta = parse_plate_metadata()
    annotation_level = 'annotation'

    for tissue in args.tissues:
        print(tissue)

        print('Load dataset')
        ds_tissue = get_dataset(
                tissue,
                membrane_only=False,
                regenerate=args.regenerate)

        if args.subtissue == 'all':
            subtissues = np.unique(ds_tissue.samplesheet['subtissue'])
        else:
            subtissues = [args.subtissue]

        for subtissue in subtissues:
            if subtissue is None:
                ds = ds_tissue
            else:
                ds = ds_tissue.query_samples_by_metadata(
                        'subtissue == @subtissue',
                        local_dict=locals(),
                        inplace=False)

            if args.cell_types is None:
                cell_types = np.unique(ds.samplesheet[annotation_level])
            else:
                cell_types = args.cell_types

            if len(args.markers) == 1:
                xname = args.markers[0]
                xgene = channels_to_genes.get(xname, xname)

                data = ds.counts.loc[[xgene]].T
                data['cell type'] = ds.samplesheet['annotation']

                fig, ax = plt.subplots()
                sns.violinplot(
                        data=data,
                        x='cell type',
                        y=xgene,
                        ax=ax,
                        zorder=5,
                        scale='width'
                        )
                ax.grid(True)
                for tk in ax.xaxis.get_ticklabels():
                    tk.set_rotation(90)
                if subtissue is not None:
                    ax.set_title('{:}, {:}'.format(tissue, subtissue))
                else:
                    ax.set_title('{:}, {:}'.format(tissue, 'all subtissues'))
                plt.tight_layout()

            if len(args.markers) == 2:
                xname, yname = args.markers
                xgene = channels_to_genes.get(xname, xname)
                ygene = channels_to_genes.get(yname, yname)

                fig, axs = make_subplots(nplots=len(cell_types))
                for ax, cell_type in zip(axs, cell_types):
                    # Set identity one VS all
                    col = 'cell_type: {:}'.format(cell_type)
                    ds.samplesheet[col] = False
                    ds.samplesheet.loc[
                            ds.samplesheet[annotation_level] == cell_type,
                            'cell_type: {:}'.format(cell_type)] = True

                    data = ds.counts.loc[[xgene, ygene]].T
                    data['cell type'] = ds.samplesheet[col]
                    c = ['steelblue' if x else 'grey' for x in data['cell type']]
                    zorder = np.where(data['cell type'].values, 5., 4.)
                    data.plot(
                            kind='scatter',
                            x=xgene,
                            y=ygene,
                            color=c,
                            s=20,
                            alpha=0.5,
                            ax=ax,
                            #FIXME
                            #zorder=zorder,
                            )
                    ax.set_title(cell_type)

    plt.ion()
    plt.show()
