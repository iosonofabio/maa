# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/10/18
content:    Test singlet config by loading the data
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import seaborn as sns



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

    hpa = parse_hpa()
