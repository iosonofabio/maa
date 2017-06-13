# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/06/17
content:    Test fcs parsing
'''
# Modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from maa.filenames import get_fcs_filenames
from maa.facs import FacsSample, FacsSort


# Script
if __name__ == '__main__':

    mouse = '30_1_M'
    tissue = 'spleen'

    #fns = get_fcs_filenames(mouse, tissue)

    ## Use tmp file for now (check the parser)
    #fn = '../data/tmp/Spleen sort PI stained - 1_MAA000032.fcs'
    #import fcsparser
    #metadata, data = fcsparser.parse(fn, meta_data_only=False, reformat_meta=True)

    s = FacsSort('MAA000032')
    data = s.get_fcs_data()
    s.plot_fcs_data(
            axes=(
                ('FSC-A', 'SSC-A'),
                ('FSC-H', 'FSC-W'),
                ('PI-A', 'FSC-A')),
            scales=(
                ('linear', 'linear'),
                ('linear', 'linear'),
                ('log', 'linear')),
            include_index_sort=True,
            )

    plt.ion()
    plt.show()

    #s = FacsSample(mouse, tissue)
    #so = s.sorts[0]

