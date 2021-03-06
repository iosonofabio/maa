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


def plot_svm(X, y, clf, ax=None):
    '''Plot SVM classifier'''
    if ax is None:
        fig, ax = plt.subplots()

    colors = plt.cm.Paired([0.0, 1.0])
    colors[0, -1] = 0.4
    colors[1, -1] = 0.7
    c = np.zeros((len(y), 4))
    c[y == 0] = colors[0]
    c[y == 1] = colors[1]

    if X.shape[1] == 1:
        from scipy.optimize import minimize_scalar
        def fun(x, offset):
            return (clf.decision_function([[x]])[0] - offset)**2
        discr = minimize_scalar(fun, args=(0,), bounds=[0, 6]).x
        dis_low = minimize_scalar(fun, args=(-0.5,), bounds=[0, 6]).x
        dis_high = minimize_scalar(fun, args=(+0.5,), bounds=[0, 6]).x
        df = pd.DataFrame([X[:, 0], y], index=['x', 'identity']).T
        sns.swarmplot(
                x='x', y='identity', data=df, ax=ax,
                orient='h',
                alpha=0.7,
                )
        ax.axvline(discr)
        ax.axvline(dis_low, ls='--')
        ax.axvline(dis_high, ls='--')

    else:
        ax.scatter(
                X[:, 0], X[:, 1],
                color=c,
                zorder=10,
                s=20)

        x_min = X[:, 0].min()
        x_max = X[:, 0].max()
        y_min = X[:, 1].min()
        y_max = X[:, 1].max()

        XX, YY = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
        Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])

        # Put the result into a color plot
        Z = Z.reshape(XX.shape)
        ax.pcolormesh(XX, YY, Z > 0, cmap=plt.cm.Paired, alpha=0.05)
        ax.contour(XX, YY, Z, colors=['k', 'k', 'k'],
                   linestyles=['--', '-', '--'], levels=[-.5, 0, .5])

    return ax


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
    parser.add_argument('--save', action='store_true',
                        help='Store to file instead of showing')
    parser.add_argument('--only-commercial', action='store_true',
                        help='Select only commercially available markers')
    parser.add_argument('--cell-types', nargs='+', default=None,
                        help='Limit to some cell types')
    parser.add_argument('--plot-suboptimal', type=int, default=0,
                        help='Plot N suboptimal clusters as well')
    parser.add_argument('--max-candidates', type=int, default=20,
                        help='Max number of candidates to try and classify as pairs')
    parser.add_argument('--go-contains', default=None,
                        help='Select genes with a substring in their GO annotation')
    parser.add_argument('--go-exclude', default=None,
                        help='Select genes missing a substring in their GO annotation')
    parser.add_argument('--go-annotate', default=None,
                        help='Annotate genes that contain this GO substring with a circle (○)')
    parser.add_argument('--exclude-genes', nargs='+', default=(),
                        help='Exclude these genes from the candidate lists')
    parser.add_argument('--kernel', default='linear',
                        choices=('linear', 'rbf', 'poly2', 'poly3', 'sigmoid'),
                        help='Kernel for the SVM classifier')
    args = parser.parse_args()

    if (len(args.tissues) == 1) and (args.tissues[0] == 'all'):
        args.tissues = tuple(tissues_prediction.keys())

    if (args.go_contains is not None) and (args.go_exclude is not None):
        raise ValueError('You can use either go-contains xor go-exclude, not both')

    if args.go_annotate is not None:
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

    for tissue in args.tissues:
        print(tissue)

        fn_cache = '../../data/cache/{:}_cache.tsv'

        print('Load dataset')
        ds = get_dataset(
                tissue,
                membrane_only=True,
                regenerate=args.regenerate,
                go_contains=args.go_contains,
                go_exclude=args.go_exclude)

        # Usually we want only a subtissue
        #FIXME: generalize this
        subtissues = np.unique(ds.samplesheet['subtissue'])
        subtissue = subtissues[0]
        ds.query_samples_by_metadata(
                'subtissue == @subtissue',
                local_dict=locals(),
                inplace=True)

        annotation_level = 'annotation'
        classifiers = []
        if args.cell_types is None:
            cell_types = np.unique(ds.samplesheet[annotation_level])
        else:
            cell_types = args.cell_types
        for cell_type in cell_types:
            # Set identity one VS all
            col = 'cell_type: {:}'.format(cell_type)
            ds.samplesheet[col] = False

            # NOTE: colon has two types of stem cells, but they are really the same
            if cell_type == 'Undiff. Cell':
                ind = ds.samplesheet[annotation_level].isin(
                        ['Cycling Undiff. Cell',
                         'Non-Cycling Undiff. Cell'])
                ds.samplesheet.loc[ind, annotation_level] = 'Undiff. Cell'

            ds.samplesheet.loc[
                    ds.samplesheet[annotation_level] == cell_type,
                    'cell_type: {:}'.format(cell_type)] = True

            # Get top membrane genes that separate the two
            dst = ds.split(col)
            res = dst[True].compare(dst[False]).rename(columns={'P-value': 'pval'})
            res['median_expression_True'] = dst[True].counts.median(axis=1)
            res['median_expression_False'] = dst[False].counts.median(axis=1)
            pval_max = 1e-3
            if args.only_commercial:
                pval_max *= 10
            candidates = res.query('pval <= @pval_max').nsmallest(columns=['pval'], n=100)

            if args.only_commercial:
                ind = np.intersect1d(candidates.index, ab_comm)
                candidates = candidates.loc[ind]

            if args.exclude_genes:
                candidates = candidates.loc[~candidates.index.isin(args.exclude_genes)]

            # Try out SVMs
            from sklearn import svm
            kernel = args.kernel
            if 'poly' in kernel:
                degree = int(kernel[-1])
                kernel = 'poly'
            else:
                degree = 3

            if len(candidates) == 0:
                print('No discriminatory genes in {:}'.format(cell_type))
                classifiers.append({
                    'classifier': None,
                    'cell type': cell_type,
                    'tissue': tissue,
                    })
                continue

            if len(candidates) == 1:
                g1 = candidates.index[0]
                print('Only one discriminatory gene in {:}: {:}'.format(cell_type, g1))
                # X is [n samples, n features]
                X = ds.counts.loc[[g1]].values.T
                y = ds.samplesheet[col].values.astype(int)
                clf = svm.SVC(C=1, kernel=kernel, degree=degree, class_weight={1: 10})
                clf.fit(X, y)

                # Check enrichment
                true_pos = ((clf.decision_function(X) > 0) & y).sum()
                false_neg = ((clf.decision_function(X) <= 0) & y).sum()
                false_pos = ((clf.decision_function(X) > 0) & (~y)).sum()
                true_neg = ((clf.decision_function(X) <= 0) & (~y)).sum()
                precision = 1.0 * true_pos / (true_pos + false_pos)
                prevalence = y.mean()
                enrichment = precision / prevalence
                # sensitivity aka recall
                recall = 1.0 * true_pos / y.sum()
                specificity = 1.0 * true_neg / (~y).sum()

                clas_best = {
                    'X': X,
                    'y': y,
                    'genes': (g1,),
                    'classifier': clf,
                    'true_pos': true_pos,
                    'false_pos': false_pos,
                    'false_neg': false_neg,
                    'true_neg': true_neg,
                    'precision': precision,
                    'prevalence': prevalence,
                    'enrichment': enrichment,
                    'recall': recall,
                    'specificity': specificity,
                    'cell type': cell_type,
                    'tissue': tissue,
                    }

            else:
                print('Classifying {:}'.format(cell_type))
                classifiers_sub = []
                for i in range(min(len(candidates), args.max_candidates)):
                    for j in range(i):
                        g1 = candidates.index[i]
                        g2 = candidates.index[j]

                        ij = len(classifiers_sub)
                        if ((ij == 0) or not (ij % 100)):
                            print('Pair n {:}'.format(ij + 1))

                        # X is [n samples, n features]
                        X = ds.counts.loc[[g1, g2]].values.T
                        y = ds.samplesheet[col].values.astype(int)
                        clf = svm.SVC(C=1, kernel=kernel, degree=degree, class_weight={1: 10})
                        clf.fit(X, y)

                        # Check enrichment
                        true_pos = ((clf.decision_function(X) > 0) & y).sum()
                        false_neg = ((clf.decision_function(X) <= 0) & y).sum()
                        false_pos = ((clf.decision_function(X) > 0) & (~y)).sum()
                        true_neg = ((clf.decision_function(X) <= 0) & (~y)).sum()
                        precision = 1.0 * true_pos / (true_pos + false_pos)
                        prevalence = y.mean()
                        enrichment = precision / prevalence
                        # sensitivity aka recall
                        recall = 1.0 * true_pos / y.sum()
                        specificity = 1.0 * true_neg / (~y).sum()

                        classifiers_sub.append({
                            'X': X,
                            'y': y,
                            'cellnames': ds.counts.columns.tolist(),
                            'genes': (g1, g2),
                            'classifier': clf,
                            'true_pos': true_pos,
                            'false_pos': false_pos,
                            'false_neg': false_neg,
                            'true_neg': true_neg,
                            'precision': precision,
                            'prevalence': prevalence,
                            'enrichment': enrichment,
                            'recall': recall,
                            'specificity': specificity,
                            'cell type': cell_type,
                            'tissue': tissue,
                            'precision+recall': precision + recall,
                            })

                from operator import itemgetter
                clas_best = max(
                        classifiers_sub,
                        key=itemgetter('precision+recall'))

                if args.plot_suboptimal > 0:
                    print('Plotting suboptimal classifiers')
                    nplots = min(args.plot_suboptimal + 1, len(classifiers_sub))
                    fig, axs = make_subplots(nplots)
                    axs = axs.ravel()
                    classifiers_sub_sorted = sorted(classifiers_sub, key=itemgetter('precision+recall'), reverse=True)
                    for (d, ax) in zip(classifiers_sub_sorted, axs):
                        g1, g2 = d['genes']
                        clf = d['classifier']
                        X = d['X']
                        y = d['y']
                        plot_svm(X, y, clf, ax=ax)
                        xlabel = 'log10 expression of {:}'.format(g1)
                        if g1 in ab_comm:
                            xlabel += '*'
                        if args.go_annotate is not None and args.go_annotate in go.loc[g1, 'GONames']:
                            xlabel += '○'
                        ax.set_xlabel(xlabel)
                        ylabel = 'log10 expression of {:}'.format(g2)
                        if g2 in ab_comm:
                            ylabel += '*'
                        if args.go_annotate is not None and args.go_annotate in go.loc[g2, 'GONames']:
                            ylabel += '○'
                        ax.set_xlabel(xlabel)
                        ax.set_ylabel(ylabel)
                        ax.grid(False)
                        ax.set_title(
                            'p={:.0%}, {:.1f}x, r={:.0%}'.format(d['precision'], d['enrichment'], d['recall']),
                            fontsize=9)

                    fig.suptitle('{:s}: prevalence={:.0%}'.format(cell_type, d['prevalence']))
                    plt.tight_layout(rect=(0, 0, 1, 0.96))

                    genes_print = set()
                    for d in classifiers_sub_sorted[:nplots]:
                        g1, g2 = d['genes']
                        genes_print.add(g1)
                        genes_print.add(g2)
                    genes_print = list(genes_print)

                    # Annotate some candidates with a star
                    if args.go_annotate is not None:
                        for i, g in enumerate(genes_print):
                            if args.go_annotate in go.loc[g, 'GONames']:
                                genes_print[i] = g+' ○'

                    print('Antibodies to test:')
                    print('\n'.join(genes_print))

                    print('Antibody pairs:')
                    for d in classifiers_sub_sorted[:nplots]:
                        g1, g2 = d['genes']
                        print('{:}\t{:}'.format(g1, g2))

            classifiers.append(clas_best)

        print('Plotting')
        nplots = len(classifiers)
        fig, axs = make_subplots(nplots)
        if len(axs) > len(classifiers):
            for ax in axs[len(classifiers):]:
                ax.axis('off')

        for (d, ax) in zip(classifiers, axs):
            clf = d['classifier']
            if clf is None:
                ax.set_title(d['cell type'], fontsize=9)
                ax.axis('off')
                continue
            if len(d['genes']) == 2:
                g1, g2 = d['genes']
            else:
                g1 = d['genes'][0]
                g2 = None
            X = d['X']
            y = d['y']
            plot_svm(X, y, clf, ax=ax)
            xlabel = 'log10 expression of {:}'.format(g1)
            if g1 in ab_comm:
                xlabel += '*'
            if args.go_annotate is not None and args.go_annotate in go.loc[g1, 'GONames']:
                xlabel += '○'
            ax.set_xlabel(xlabel)
            if g2 is not None:
                ylabel = 'log10 expression of {:}'.format(g2)
                if g2 in ab_comm:
                    ylabel += '*'
                if args.go_annotate is not None and args.go_annotate in go.loc[g2, 'GONames']:
                    xlabel += '○'
            else:
                ylabel = ''
            ax.set_ylabel(ylabel)
            ax.grid(False)
            ax.set_title(
                    '{:s}: p={:.0%}→{:.0%} ({:.1f}x), r={:.0%}'.format(
                        d['cell type'], d['prevalence'], d['precision'], d['enrichment'], d['recall']),
                    fontsize=9)
        fig.suptitle(tissue)
        plt.tight_layout(rect=(0, 0, 1, 0.95))

        if args.save:
            from sklearn.externals import joblib
            import tarfile
            import json
            fields = (
                    'tissue',
                    'precision',
                    'recall',
                    'enrichment',
                    'prevalence',
                    'specificity',
                    'cell type',
                    'genes',
                    )

            for classifier in classifiers:
                if subtissue is not None:
                    fn_glb = '../../data/classifiers/{:}_{:}_{:}_antibodies'.format(
                            tissue.lower(),
                            subtissue.lower(),
                            classifier['cell type'].replace(' ', '_'),
                            )
                else:
                    fn_glb = '../../data/classifiers/{:}_{:}_antibodies'.format(
                            tissue.lower(),
                            classifier['cell type'].replace(' ', '_'),
                            )
                fn_model = fn_glb+'.model.pickle'
                fn_train = fn_glb+'.train.npz'
                fn_meta = fn_glb+'.metadata.json'
                fn_bundle = fn_glb+'.tar.gz'

                # Save classifier
                clf = classifier['classifier']
                joblib.dump(clf, fn_model)

                # Save metadata
                meta = {k: classifier[k] for k in fields}
                with open(fn_meta, 'wt') as f:
                    json.dump(meta, f)

                # Save training
                np.savez_compressed(
                        fn_train,
                        X=classifier['X'],
                        y=classifier['y'],
                        cellnames=classifier['cellnames'])

                # Bundle up
                with tarfile.open(fn_bundle, 'w:gz') as f:
                    f.add(fn_model, arcname=os.path.basename(fn_model))
                    f.add(fn_meta, arcname=os.path.basename(fn_meta))
                    f.add(fn_train, arcname=os.path.basename(fn_train))

        plt.ion()
        plt.show()
