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


if __name__ == '__main__':

    pa = argparse.ArgumentParser(
        description='Test loading various tabula muris tissues',
        )
    pa.add_argument(
        '--tissues', default=tissues,
        nargs='+',
        choices=tissues,
        help='Tissue of origin',
        )
    args = pa.parse_args()

    dss = {}
    for tissue in args.tissues:
        print(tissue)
        ds = Dataset(
            counts_table=tissue,
            samplesheet='alltissues',
            featuresheet='alltissues',
            )
        dss[tissue] = ds
