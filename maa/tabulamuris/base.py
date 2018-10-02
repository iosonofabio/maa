# vim: fdm=indent
'''
author:     Fabio Zanini
date:       02/10/18
content:    Basic tools for working with the Tabula Muris data
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr

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


def load_hits(tissue, method='KS', cell_type=None):
    '''Load the table of hits for a tissue and possibly cell type

    Args:
        tissue (str): the tissue to analyze
        method (str): the method used to find hits
        cell_type (str): restrict the list to a specific cell type

    Returns:
        hits (pd.DataFrame): the table of hits.
    '''
    fn = '../../data/tabula_muris/FACS/hits/hits_{:}_{:}.tsv'.format(
            method,
            tissue,
            )
    hits = pd.read_csv(fn, sep='\t', index_col=0)
    if cell_type is not None:
        hits.query('cell_type == @cell_type', inplace=True)
    return hits


def load_mouse_human_homology():
    # MIA2 in humans and Try5 in mouse are present twice:
    # - MIA2 one of the entries is very incomplete, cutting it
    # - Try5 both entries are exactly identical
    # I deleted the 2nd entry for both from the file
    fn = '../../data/mouse_human_homology/mouse_human.tsv'
    mho = pd.read_csv(fn, sep='\t', comment='#')
    mho['Organism'] = [x.split(',')[0] for x in mho['Common Organism Name']]
    mho.set_index(['Symbol', 'Organism'], drop=False, inplace=True)
    mho.rename(
        columns={x: x.replace(' ', '') for x in mho.columns},
        inplace=True)
    return mho


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
