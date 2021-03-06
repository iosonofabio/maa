# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/11/17
content:    Merge scattering and antibody-based predictors.
'''
# Modules
import os
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.externals import joblib
import json

import matplotlib.pyplot as plt
import seaborn as sns


# Functions
class CombinedClassifier:
    def __init__(self, classifiers, logic='AND'):
        self.classifiers = classifiers
        self.logic = logic

    def predict(self, Xs):
        y = (self.classifiers[0].predict(Xs[0]) > 0)
        if len(self.classifiers) > 1:
            if self.logic == 'AND':
                for clf, X in zip(self.classifiers[1:], Xs[1:]):
                    y &= (clf.predict(X) > 0)
            elif self.logic == 'OR':
                for clf, X in zip(self.classifiers[1:], Xs[1:]):
                    y |= (clf.predict(X) > 0)
            else:
                raise ValueError('Combination logic not understood: {}'.format(self.logic))

        return (y * 2) - 1


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


def plot_classifier_antibody(X, y, clf, ax=None):
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


def plot_classifier_scattering(X, y, clf, ax=None):
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



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('tissues', nargs='+',
                        help='tissues to study')
    parser.add_argument('--cell-types', nargs='+', required=True,
                        help='Limit to some cell types')
    parser.add_argument('--subtissue', default=None,
                        help='Limit to a subtissue. To split by subtissue use "all"')
    parser.add_argument('--save', action='store_true',
                        help='Store to file instead of showing')
    parser.add_argument('--save-website', action='store_true',
                        help='Save result for the website JSON)')
    parser.add_argument('--combination-logic', default='AND',
                        choices=['AND', 'OR'],
                        help='Combination logic between scattering and antibodies')
    args = parser.parse_args()

    # Get the list of commercially available antibodies
    ab_comm = []
    # Biolegend
    ab_comm_table = parse_biolegend()
    ab_unique = np.unique(ab_comm_table.dropna(subset=['GeneName'], axis=0)['GeneName'])
    ab_comm.append(ab_unique)
    # TODO: other vendors
    if len(ab_comm):
        ab_comm = np.unique(np.concatenate(ab_comm))


    for tissue in args.tissues:
        classifiers = []
        for cell_type in args.cell_types:
            print(cell_type)
            clfs = {}
            for clf_type in ('scattering', 'antibodies', 'combined'):
                if args.subtissue is not None:
                    fn_glb = '../../data/classifiers/{:}_{:}_{:}_{:}'.format(
                            tissue.lower(),
                            args.subtissue.lower(),
                            cell_type.replace(' ', '_'),
                            clf_type,
                            )
                else:
                    fn_glb = '../../data/classifiers/{:}_{:}_{:}'.format(
                            tissue.lower(),
                            cell_type.replace(' ', '_'),
                            clf_type,
                            )

                fn_model = fn_glb+'.model.pickle'
                fn_train = fn_glb+'.train.npz'
                fn_meta = fn_glb+'.metadata.json'
                fn_bundle = fn_glb+'.tar.gz'

                with open(fn_meta, 'r') as f:
                    classifier = json.load(f)
                clf = joblib.load(fn_model)
                classifier['classifier'] = clf
                train = np.load(fn_train)
                classifier['X'] = train['X']
                classifier['y'] = train['y']
                classifier['cellnames'] = train['cellnames']

                clfs[clf_type] = classifier
            classifiers.append(clfs)

            # Combine the classifiers
            cells_common = np.intersect1d(
                    clfs['scattering']['cellnames'],
                    clfs['antibodies']['cellnames'],
                    )
            Xs = pd.DataFrame(
                    data=clfs['scattering']['X'],
                    index=clfs['scattering']['cellnames'],
                    ).loc[cells_common].values
            ys = pd.Series(
                    data=clfs['scattering']['y'],
                    index=clfs['scattering']['cellnames'],
                    ).loc[cells_common].values
            clas = clfs['scattering']['classifier']
            Xa = pd.DataFrame(
                    data=clfs['antibodies']['X'],
                    index=clfs['antibodies']['cellnames'],
                    ).loc[cells_common].values
            ya = pd.Series(
                    data=clfs['antibodies']['y'],
                    index=clfs['antibodies']['cellnames'],
                    ).loc[cells_common].values
            claa = clfs['antibodies']['classifier']

            if (ys != ya).any():
                raise ValueError('The true cell identity should be the same!')

            cus = CombinedClassifier(
                    classifiers=[clas, claa],
                    logic=args.combination_logic)


        from sklearn.externals import joblib
        import json
        for clfs in classifiers:
            clf = clfs['combined']
            fn = '{:}/university/postdoc/facsweb/app/facsweb/shared/static/data/merged_predictor_{:}_{:}.json'.format(
                    os.getenv('HOME'),
                    clf['tissue'],
                    clf['cell type'].replace(' ', '_'),
                    )

            d = {}
            d['tissue'] = clf['tissue']
            d['cell type'] = clf['cell type']

            d['data'] = {}
            d['data']['scattering'] = clf['Xs'].tolist()
            d['data']['antibodies'] = clf['Xa'].tolist()
            d['data']['identity'] = clf['y'].tolist()
            d['data']['cellnames'] = clf['cellnames'].tolist()
            d['data']['scattering_axis_labels'] = [clf['xlabel'], clf['ylabel']]
            d['data']['antibody_axis_labels'] = clf['genes']
            d['data']['xlim_scattering'] = [0, clf['Xs'][:, 0].max() * 1.]
            d['data']['ylim_scattering'] = [
                    np.floor(clf['Xs'][:, 1].min()),
                    np.ceil(clf['Xs'][:, 1].max()),
                    ]
            d['data']['xlim_antibodies'] = [
                    np.floor(clf['Xa'][:, 0].min()),
                    np.ceil(clf['Xa'][:, 0].max()),
                    ]
            d['data']['ylim_antibodies'] = [
                    np.floor(clf['Xa'][:, 1].min()),
                    np.ceil(clf['Xa'][:, 1].max()),
                    ]

            d['models'] = {
                'combined': {},
                'scattering': {},
                'antibodies': {},
                }
            d['models']['combined']['precision'] = clf['precision']
            d['models']['combined']['recall'] = clf['recall']
            d['models']['combined']['logic'] = clf['combination_logic']

            # Find roots of the classifiers
            for clfname in ('scattering', 'antibodies'):
                clfi = clfs[clfname]
                xlim = [clf['X'+clfname[0]][:, 0].min(), clf['X'+clfname[0]][:, 0].max()]
                ylim = [clf['X'+clfname[0]][:, 1].min(), clf['X'+clfname[0]][:, 1].max()]
                xx = np.linspace(xlim[0], xlim[1], 500)
                yy = np.linspace(ylim[0], ylim[1], 500)
                xv, yv = np.meshgrid(xx, yy)
                grid = np.vstack([xv.ravel(), yv.ravel()]).T
                dec = clfi['classifier'].decision_function(grid)

                roots = grid[np.abs(dec) < 0.02]
                d['models'][clfname]['roots'] = roots.tolist()

                roots_pos = grid[np.abs(dec - 0.25) < 0.02]
                d['models'][clfname]['roots_pos'] = roots_pos.tolist()

                roots_neg = grid[np.abs(dec + 0.25) < 0.02]
                d['models'][clfname]['roots_neg'] = roots_neg.tolist()

                d['models'][clfname]['precision'] = clfi['precision']
                d['models'][clfname]['recall'] = clfi['recall']

            with open(fn, 'w') as f:
                json.dump(d, f)

            print('Saved to file: {:}'.format(fn))

        plt.ion()
        plt.show()
