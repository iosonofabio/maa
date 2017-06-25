# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/06/17
content:    Try to reverse engineer the Sony expdat format.
'''
# Modules
import os
import sys
import fnmatch
import re
import zipfile
import io
import argparse
import numpy as np
import pandas as pd

from maa.sony_database.sony_database import SonyExpdat



# Script
if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Get metadata out of Sony expdata.')
    paa = pa.add_argument
    paa('mice', nargs='+',
        help='Mice to export')
    paa('--output',
        help='Output to file')

    args = pa.parse_args()

    data = []
    for mice_day in args.mice:
        e = SonyExpdat(mice_day)
        e.get_filename()
        datum = e.get_index_metadata()
        data.append(datum)
    data = pd.concat(data)

    data_export = data.loc[:, [
            'Well', 'SortGate', 'Mouse', 'Tissue', 'Subtissue', 'PlateBarcode',
            'ExpId', 'Nozzle', 'CytometerLocation', 'Tube_Name', 'SortType',
            ]]

    if args.output:
        data.to_csv(args.output, sep='\t', header=True, index=False)
