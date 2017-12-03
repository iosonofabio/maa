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
from operator import itemgetter
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


def parse_facs_plate(tissue, plate, with_fcs=True):
    sorter = config[tissue].get('sorter', 'Sony')
    glb = config[tissue].get('facs_glob', '*')
    if sorter == 'Sony':
        return parse_facs_plate_sony(plate, glb=glb, with_fcs=with_fcs)
    else:
        return parse_facs_plate_aria(plate, with_fcs=with_fcs)


def parse_facs_plate_aria(plate, glb='*', with_fcs=True):
    import glob
    import fcsparser

    out = {}

    fdn = '../../data/MACAFACS/index_fcs_{:}/'.format(glb)

    if with_fcs:
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

    if with_fcs:
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


def parse_facs_plate_sony(plate, glb='*', with_fcs=True):
    import glob
    import fcsparser

    out = {}

    fdn = '../../data/MACAFACS/index_fcs_{:}/'.format(glb)

    if with_fcs:
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

    if with_fcs:
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
    if with_fcs:
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


def plot_classifier(X, y, clf, ax=None, method=None):
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
        if method == 'svm':
            def fun(x, offset):
                return (clf.decision_function([[x]])[0] - offset)**2
        elif method == 'logistic':
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
        if method == 'svm':
            Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])
        elif method == 'logistic':
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


def get_facs_data(tissue, cell_types, plates, maxplates=-1, with_fcs=True):
    facs_data = {}
    ip = 0
    for plate in plates:
        print('Parse FACS data: plate {:}'.format(plate))

        # Some plates are missing
        try:
            facs_datum = parse_facs_plate(tissue, plate, with_fcs=with_fcs)
        except IOError:
            print('Missing')
            continue

        ind_bool = facs_datum['index_data'].index.isin(cell_types.index)
        ind = facs_datum['index_data'].index[ind_bool]
        facs_datum['index_data']['Good quality'] = ind_bool

        facs_datum['index_data']["annotation"] = "missing"
        facs_datum['index_data'].loc[ind, "annotation"] = \
            cell_types.loc[ind, "annotation"]

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

        ip += 1
        if ip == maxplates:
            break

    return facs_data


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('tissues', nargs='+',
                        help='tissues to study')
    parser.add_argument('--regenerate', action='store_true',
                        help='Regenerate counts cache')
    parser.add_argument('--save', action='store_true',
                        help='Store to file instead of showing')
    parser.add_argument('--skip-fcs', action='store_true',
                        help='Skip whole FCS for faster execution')
    parser.add_argument('--maxplates', type=int, default=0,
                        help='Max number of plates to plot')
    parser.add_argument('--subtissue', default=None,
                        help='Limit to a subtissue. To split by subtissue use "all"')
    parser.add_argument('--cell-types', nargs='+', default=None,
                        help='Limit to some cell types')
    parser.add_argument('--method', default='svm',
                        choices=('svm', 'logistic'),
                        help='Method to use')
    parser.add_argument('--kernel', default='linear',
                        choices=('linear', 'rbf', 'poly2', 'poly3', 'sigmoid'),
                        help='Kernel for the SVM classifier')
    parser.add_argument('--class-weight', type=float, nargs=2, default=(1, 10),
                        help='Weights for the negative and positive class')
    parser.add_argument('--merge-plates', action='store_true',
                        help='Merge plates from the same tissue or subtissue')
    args = parser.parse_args()

    if (len(args.tissues) == 1) and (args.tissues[0] == 'all'):
        args.tissues = tuple(tissues_prediction.keys())

    plate_meta = parse_plate_metadata()

    for tissue in args.tissues:
        print(tissue)

        cell_types, plates = parse_annotations(tissue)
        facs_data = get_facs_data(
                tissue,
                cell_types,
                plates,
                maxplates=args.maxplates,
                with_fcs=not args.skip_fcs)

        annotation_level = 'annotation'
        cell_types_unique = np.unique(cell_types[annotation_level])

        if args.subtissue == 'all':
            subtissues = np.unique(plate_meta.loc[plates, 'subtissue'])
        else:
            subtissues = [args.subtissue]

        for subtissue in subtissues:
            index_data_st = []
            for plate, facs_datum in facs_data.items():
                subtissue_plate = plate_meta.loc[plate, 'subtissue']
                if (subtissue is not None) and (subtissue != subtissue_plate):
                    continue

                print('Plate {:}, {:}, {:}'.format(plate, tissue, subtissue_plate))
                facs_datum = facs_data[plate]
                index_datap = facs_datum['index_data']
                index_datap['SSC-A-log'] = 0.1 + np.log10(index_datap['SSC-A'])
                index_datap['FSC-A-1e5'] = index_datap['FSC-A'] / 1e5
                index_datap['plate'] = plate
                index_datap['subtissue'] = subtissue_plate

                index_data_st.append(index_datap)

            if args.merge_plates:
                index_data_st = [pd.concat(index_data_st, axis=0)]

            for index_data in index_data_st:
                plate = ', '.join(np.unique(index_data['plate']))
                subtissue_plate = ', '.join(np.unique(index_data['subtissue']))

                from sklearn import svm
                from sklearn import linear_model
                if args.method == 'logistic':
                    pass
                elif args.method == 'svm':
                    kernel = args.kernel
                    if 'poly' in kernel:
                        degree = int(kernel[-1])
                        kernel = 'poly'
                    else:
                        degree = 3
                else:
                    raise ValueError('Method not recognized')

                classifiers = []
                if args.cell_types is None:
                    cell_types_plate = np.unique(index_data[annotation_level])
                else:
                    cell_types_plate = args.cell_types

                for cell_type in cell_types_plate:
                    if cell_type == 'missing':
                        continue
                    print('Classifying {:}'.format(cell_type))
                    col = 'cell_type: {:}'.format(cell_type)

                    index_data.loc[:, col] = False
                    index_data.loc[
                            index_data[annotation_level] == cell_type,
                            'cell_type: {:}'.format(cell_type)] = True

                    # Cells that did not make it to the sequencer, i.e. missing,
                    # can trip up the classifier
                    index_data = index_data.loc[index_data[annotation_level] != 'missing']

                    # X is [n samples, n features]
                    X = index_data.loc[:, ['FSC-A-1e5', 'SSC-A-log']].values
                    y = index_data[col].values.astype(int)

                    if tuple(args.class_weight) != (0, 0):
                        class_weights = [args.class_weight]
                    else:
                        class_weights = [
                                (1, 3),
                                (1, 5),
                                (1, 7.5),
                                (1, 10),
                                (1, 30),
                                (1, 50),
                                (1, 100),
                                ]

                    classifiers_cw = []
                    for class_weight in class_weights:
                        if args.method == 'logistic':
                            clf = linear_model.LogisticRegression(
                                    C=1,
                                    class_weight={
                                        0: class_weight[0],
                                        1: class_weight[1]})
                            clf.fit(X, y)
                        elif args.method == 'svm':
                            clf = svm.SVC(
                                    C=1, kernel=kernel, degree=degree,
                                    class_weight={
                                        0: class_weight[0],
                                        1: class_weight[1]})
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

                        classifiers_cw.append({
                            'X': X,
                            'y': y,
                            'cellnames': index_data.index.tolist(),
                            'xlabel': 'FSC-A-1e5',
                            'ylabel': 'SSC-A-log',
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
                            'precision+recall': 2 * precision + recall,
                            })

                    classif = max(
                            classifiers_cw,
                            key=itemgetter('precision+recall'),
                            )
                    classifiers.append(classif)

                nplots = len(cell_types_plate)
                fig, axs = make_subplots(nplots)
                if len(axs) > len(classifiers):
                    for ax in axs[len(classifiers):]:
                        ax.axis('off')

                for (d, ax) in zip(classifiers, axs):
                    clf = d['classifier']
                    X = d['X']
                    y = d['y']
                    plot_classifier(X, y, clf, ax=ax, method=args.method)
                    ax.set_xlabel('FSC-1e5')
                    ax.set_ylabel('SSC-A-log')
                    ax.grid(False)
                    ax.set_title(
                            '{:s}: p={:.0%}→{:.0%} ({:.1f}x), r={:.0%}'.format(
                                d['cell type'], d['prevalence'], d['precision'], d['enrichment'], d['recall']),
                            fontsize=9)

                fig.suptitle('{:}, {:}, plate {:}'.format(tissue, subtissue_plate, plate))
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
                        'xlabel',
                        'ylabel',
                        )

                for classifier in classifiers:
                    if subtissue is not None:
                        fn_glb = '../../data/classifiers/{:}_{:}_{:}_scattering'.format(
                                tissue.lower(),
                                subtissue.lower(),
                                classifier['cell type'].replace(' ', '_'),
                                )
                    else:
                        fn_glb = '../../data/classifiers/{:}_{:}_scattering'.format(
                                tissue.lower(),
                                classifier['cell type'].replace(' ', '_'),
                                )

                    fn_model = fn_glb+'.model.pickle'
                    fn_meta = fn_glb+'.metadata.json'
                    fn_train = fn_glb+'.train.npz'
                    fn_bundle = fn_glb+'.tar.gz'

                    # Save classifier
                    clf = classifier['classifier']
                    joblib.dump(clf, fn_model)

                    # Save training
                    np.savez_compressed(
                            fn_train,
                            X=classifier['X'],
                            y=classifier['y'],
                            cellnames=classifier['cellnames'])

                    # Save metadata
                    meta = {k: classifier[k] for k in fields}
                    with open(fn_meta, 'wt') as f:
                        json.dump(meta, f)

                    # Bundle up
                    with tarfile.open(fn_bundle, 'w:gz') as f:
                        f.add(fn_model, arcname=os.path.basename(fn_model))
                        f.add(fn_train, arcname=os.path.basename(fn_train))
                        f.add(fn_meta, arcname=os.path.basename(fn_meta))

        plt.ion()
        plt.show()
