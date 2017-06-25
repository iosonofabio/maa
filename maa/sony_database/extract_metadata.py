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


# Functions / classes
class SonyForGoogleMeta(SonyExpdat):

    @staticmethod
    def convert_data_for_google_sheet(data):
        cols = ['Plate.barcode',
                'Lysis Plate Batch',
                'dNTP.batch',
                'oligodT.order.no',
                'plate.type',
                'preparation.site',
                'date.prepared',
                'date.sorted',
                'tissue',
                'subtissue',
                'Mouse ID (age_#_sex)',
                'FACS.selection',
                'nozzle.size',
                'FACS.instument',
                'Experiment ID',
                'Columns sorted',
                'Double check',
                'Plate number']

        data = data.copy()
        data['Lysis Plate Batch'] = ''
        data['dNTP.batch'] = ''
        data['oligodT.order.no'] = ''
        data['preparation.site'] = 'Biohub'
        data['date.prepared'] = ''
        data['FACS.selection'] = ''
        data['Columns sorted'] = ''
        data['Double check'] = ''
        data['Plate number'] = ''
        data.rename(columns={
            'PlateBarcode': 'Plate.barcode',
            'PlateType': 'plate.type',
            'CreateDate': 'date.sorted',
            'Tissue': 'tissue',
            'Subtissue': 'subtissue',
            'Mouse': 'Mouse ID (age_#_sex)',
            'Nozzle': 'nozzle.size',
            'CytometerLocation': 'FACS.instument',
            'ExpId': 'Experiment ID',
            }, inplace=True)
        data.set_index('Plate.barcode', inplace=True, drop=False)
        data.sort_index(inplace=True)

        # Reformat dates (clearly)
        for ix, datum in data.iterrows():
            data.at[ix, 'date.sorted'] = datum['date.sorted'].split()[0].replace('-', '')[2:]

        return data[cols]


# Script
if __name__ == '__main__':

    pa = argparse.ArgumentParser(description='Get metadata out of Sony expdata.')
    paa = pa.add_argument
    paa('mice', nargs='+',
        help='Mice to export')
    paa('--output',
        help='Output to file')

    args = pa.parse_args()

    #mice_pairs = ('1_6_M 1_7_M', '3_8_M 3_9_M', '3_10_M 3_11_M')
    data = []
    for mice_day in args.mice:
        e = SonyForGoogleMeta(mice_day)
        e.get_filename()
        datum = e.get_metadata()
        data.append(datum)
    data = pd.concat(data)
    data_gs = e.convert_data_for_google_sheet(data)

    if args.output:
        data_gs.to_csv(args.output, sep='\t', header=True, index=False)
