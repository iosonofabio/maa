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

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet import Dataset


tissues = [
    'aorta',
    'bladder',
    'brain_myeloid',
    'brain_nonmyeloid',
    'diaphragm',
    'fat',
    'heart',
    'kidney',
    'gut',
    'liver',
    'lung',
    'mammary',
    'marrow',
    'muscle',
    'pancreas',
    'skin',
    'spleen',
    'thymus',
    'tongue',
    'trachea',
    ]


class DatasetTM(Dataset):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_samples_annotated()
        self.set_coverage()

    def set_samples_annotated(self):
        is_annotated = []
        for cell, row in self.samplesheet.iterrows():
            if isinstance(row['cell_ontology_class'], str):
                is_annotated.append(True)
            else:
                is_annotated.append(False)
        self._samplesheet['is_annotated'] = is_annotated

    def set_coverage(self):
        self._samplesheet['coverage'] = self.counts.sum(axis=0)

    @staticmethod
    def parse_hpa(skip_doublets=True):
        '''Parse Human Protein Atlas TSV table

        Args:
            skip_doublets (bool): Whether to skip genes with two entries. Those are
            alternative splicings and around ~10 genes including COG8.
        '''
        from collections import Counter

        fn = '../../data/human_protein_atlas/proteinatlas.tsv'
        hpa = pd.read_csv(
            fn,
            sep='\t',
            )
        if not skip_doublets:
            return hpa
        counts = Counter(hpa['Gene'].tolist())
        fnu = [f for f, val in counts.items() if val > 1]
        hpa = hpa.loc[~hpa['Gene'].isin(fnu)]
        return hpa.set_index('Gene', drop=False)


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

        hits = (hits.iloc[:500]
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
        hits = hits.loc[['membrane' in x for x in hits['HPA_Sub']]]
        print(hits[['rank', '+/-', 'diff', 'HPA_Rel_IH', 'HPA_Rel_IF', 'HPA_Sub']])
        print()
        print()
