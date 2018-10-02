# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/10/18
content:    Test identification of surface markers by tissue and cell type.
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import seaborn as sns

from maa.tabulamuris.base import DatasetTM, tissues


def find_hits_KS(ds, save=True):
    '''Find hits for biomarkers via Kolmogorov Smirnov.

    Args:
        ds (DatasetTM): the dataset to analyze
        save (bool): Whether to save the hit list to file.

    Returns:
        hits (pd.Dataframe): the table of hits for all cell types.

    This is not expected to give the final list of biomarkers, but just
    selects from the jungle of features the ones that are somewhat likely
    to make the cut. More sophisticated downstream criteria are used after
    this script.
    '''
    cell_types = np.unique(ds.samplesheet['cell_ontology_class'])

    hits_alltypes = []
    for ct in cell_types:
        print(ct)
        dss = ds.split('is {:}'.format(ct))
        comp = dss[True].compare(dss[False])
        comp.rename(columns={'P-value': 'KS'}, inplace=True)

        dss[True].counts.log(inplace=True)
        dss[False].counts.log(inplace=True)
        mt = dss[True].counts.get_statistics(metrics=['mean'])
        mf = dss[False].counts.get_statistics(metrics=['mean'])
        comp['logmean+'] = mt['mean']
        comp['logmean-'] = mf['mean']
        comp['diff'] = comp['logmean+'] - comp['logmean-']
        comp['localization'] = ds.featuresheet['GO term name']

        # Assign KS rank
        hits = comp.sort_values('KS')
        hits['rank'] = np.arange(len(hits)) + 1
        hits['+/-'] = ['+' if x > 0 else '-' for x in hits['diff']]

        hits = (hits.iloc[:1000]
                    .sort_values('diff', ascending=False)
                    .query('localization == "membrane"'))

        # Cross-check the HPA
        hits['HPA_Gene'] = ''
        hits['HPA_Antibody'] = ''
        hits['HPA_Rel_IH'] = ''
        hits['HPA_Rel_IF'] = ''
        hits['HPA_Sub'] = ''
        for hit in hits.index:
            # Gene name conversion mouse -> human
            human_gene = hit.upper()
            if human_gene in hpa.index:
                hits.loc[hit, 'HPA_Gene'] = human_gene
                hits.loc[hit, 'HPA_Antibody'] = hpa.loc[human_gene, 'Antibody']
                hits.loc[hit, 'HPA_Rel_IH'] = hpa.loc[human_gene, 'Reliability (IH)']
                hits.loc[hit, 'HPA_Rel_IF'] = hpa.loc[human_gene, 'Reliability (IF)']
                sl = hpa.loc[human_gene, 'Subcellular location']
                if isinstance(sl, str):
                    sl = sl.replace('<br>', ', ')
                else:
                    sl = ''
                hits.loc[hit, 'HPA_Sub'] = sl

        # Only select hits that are also plasma membrane in the HPA
        hits = hits.loc[['Plasma membrane' in x for x in hits['HPA_Sub']]]
        print(hits[['rank', '+/-', 'diff', 'HPA_Rel_IH', 'HPA_Rel_IF', 'HPA_Sub']])
        print()
        print()

        hits['cell_type'] = ct
        hits['tissue'] = args.tissue

        hits_alltypes.append(hits)
    hits_alltypes = pd.concat(hits_alltypes, axis=0)
    if save:
        fn = '../../data/tabula_muris/FACS/hits/hits_KS_{:}.tsv'.format(
                args.tissue)
        hits_alltypes.to_csv(fn, sep='\t', index=True)

    return hits_alltypes


if __name__ == '__main__':

    pa = argparse.ArgumentParser(
        description='Test loading various tabula muris tissues',
        )
    pa.add_argument(
        '--tissue',
        required=True,
        choices=tissues,
        help='Tissue of origin',
        )
    args = pa.parse_args()

    print('Load data for {:}'.format(args.tissue))
    ds = DatasetTM(
        counts_table=args.tissue,
        samplesheet='alltissues',
        featuresheet='alltissues',
        )

    print('Parse Human Protein Atlas')
    hpa = ds.parse_hpa()

    print('QC cells')
    ds.query_samples_by_metadata(
        '(is_annotated == True) and (coverage >= 50000)',
        inplace=True)

    print('Normalize')
    ds.counts.normalize(inplace=True)

    print('Split by cell type')
    cell_types = np.unique(ds.samplesheet['cell_ontology_class'])
    for ct in cell_types:
        ds.samplesheet['is {:}'.format(ct)] = ds.samplesheet['cell_ontology_class'] == ct

    print('Find hits by Kolmogorov Smirnov')
    hits_alltypes = find_hits_KS(ds)
