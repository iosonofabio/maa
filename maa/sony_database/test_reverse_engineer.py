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


# Functions / classes
class SonyExpdat:

    def __init__(self, mice, cytometer='Clark'):
        self.mice = mice
        self.cytometer = cytometer

    def get_filename(self):
        import glob
        fdn_expdat = '../data/cytometer files/'
        return glob.glob(fdn_expdat+'*'+self.mice+'/'+'*'+self.mice+'*'+self.cytometer+'.expdat')[0]

    def get_full_metadata(self):
        with zipfile.ZipFile(self.get_filename(), 'r') as zf:
            self._zf = zf

            # Let's get the experiments first (top down)
            exps = self.get_metadata(
                    'Experiments/*.ex',
                    fields=('Id', 'Name', 'Investigator', 'CreateDate'))

            # Then we get the sample groups
            sample_groups = self.get_metadata(
                    'SampleGroups/*.sg',
                    fields=('Id', 'Name', 'CreateDate', 'Experiment_Id'))

            # Then we get the tubes
            tubes = self.get_metadata(
                    'Tubes/*.t',
                    fields=('Id', 'Name', 'CreateDate', 'SampleGroup_Id'))


            # Then, we get the data sources that could have the plate barcode
            data_sources = self.get_metadata(
                    'DataSources/*.ds',
                    fields=('Id', 'Name', 'CreateDate', 'SortResult_Id', 'Tube_Id'))

            del self._zf

        # Ok, now we have all the info, we just make a flat table
        data_sources['Tube_Name'] = ''
        data_sources['SampleGroup_Name'] = ''
        data_sources['SampleGroup_Id'] = ''
        data_sources['Experiment_Name'] = ''
        data_sources['Experiment_Id'] = ''
        for ix, data_source in data_sources.iterrows():
            tube = tubes.loc[data_source['Tube_Id']]
            sample_group = sample_groups.loc[tube['SampleGroup_Id']]
            experiment = exps.loc[sample_group['Experiment_Id']]

            data_sources.at[ix, 'Tube_Name'] = tube['Name']
            data_sources.at[ix, 'SampleGroup_Name'] = sample_group['Name']
            data_sources.at[ix, 'SampleGroup_Id'] = sample_group.name
            data_sources.at[ix, 'Experiment_Name'] = experiment['Name']
            data_sources.at[ix, 'Experiment_Id'] = experiment.name

        # Add annotations
        data_sources['PlateBarcode'] = ''
        data_sources['Mouse'] = ''
        data_sources['Tissue'] = ''
        data_sources['Subtissue'] = ''
        data_sources['PlateType'] = ''
        data_sources['ExpId'] = ''
        for ix, data_source in data_sources.iterrows():
            exname = data_source['Experiment_Name']
            if exname.startswith('MAA '):
                exname = exname[4:]

            # Find plate barcode and type
            platename = data_source['Name'].replace(' ', '').rstrip('\\')
            if 'MAA' in platename:
                platetype = 'Biorad HSP3841'
            elif platename.startswith('D0'):
                platetype = 'Biorad HSP3905'
            else:
                platetype = ''
                #raise ValueError('Plate type not found for plate: '+str(data_source))

            # Find tissue name
            if 'granulocyte' in exname.lower():
                tissue = 'bone marrow'
                subtissue = 'granulocytes'
            elif 'Bcell' in exname.replace(' ', '') or 'BCell' in exname.replace(' ', ''):
                tissue = 'bone marrow'
                subtissue = 'B cells'
            elif 'Tcell' in exname.replace(' ', '') or 'TCell' in exname.replace(' ', ''):
                tissue = 'bone marrow'
                subtissue = 'T cells'
            elif 'KLS'in exname:
                tissue = 'bone marrow'
                subtissue = 'KLS'
            else:
                tissue = exname.split()[-1].lower()
                if tissue == 'colon':
                    if 'prox' in data_source['Tube_Name'].lower():
                        subtissue = 'proximal'
                    elif 'dist' in data_source['Tube_Name'].lower():
                        subtissue = 'distal'
                    else:
                        raise ValueError('Subtissue of colon not found')
                elif tissue == 'pancreas':
                    if 'endo' in data_source['Tube_Name'].lower():
                        subtissue = 'endocrine'
                    elif 'exo' in data_source['Tube_Name'].lower():
                        subtissue = 'exocrine'
                    else:
                        raise ValueError('Subtissue of pancreas not found')
                else:
                    subtissue = 'NA'

            # Find out mouse name
            exname_mice = re.findall('\d+_\d+_[MF]', exname)
            tubename_mice = re.findall('\d+_\d+_[MF]', data_source['Tube_Name'])
            if len(exname_mice) == 1:
                mouse_name = exname_mice[0]
            elif len(tubename_mice) == 1:
                mouse_name = tubename_mice[0]
            elif 'Test' in data_source['Tube_Name']:
                mouse_name = 'None_0'
            else:
                raise ValueError('Mouse name not found for this plate: '+str(data_source))

            data_sources.at[ix, 'PlateBarcode'] = platename
            data_sources.at[ix, 'Mouse'] = mouse_name
            data_sources.at[ix, 'Tissue'] = tissue
            data_sources.at[ix, 'Subtissue'] = subtissue
            data_sources.at[ix, 'PlateType'] = platetype
            data_sources.at[ix, 'ExpId'] = 'exp'+str(1+int(mouse_name.split('_')[1]))
            #data_sources.at[ix, 'BarcodeFull'] = platename+' '+mouse_name+' '+tissue+' '+subtissue

        # These are fixed at Clark
        # FIXME
        data_sources['Nozzle'] = '100'
        data_sources['CytometerLocation'] = self.cytometer

        # Filter out test sorts (not actual plates)
        data_sources = data_sources.loc[data_sources['Mouse'] != 'None_0']

        return data_sources

    def get_columns_format(self, extension, fields=None):
        '''Get meaning of the columns from hte fmt files'''
        fn = 'Formats/'+extension+'.fmt'
        with self._zf.open(fn, 'r') as f:
            lines = f.read().decode().split('\n')
        d = {}
        for line in lines:
            if 'SQLCHAR' not in line:
                continue
            line = line.rstrip('\r').split()
            col = int(line[0]) - 1
            name = line[6]
            if (fields is None) or (name in fields):
                d[col]= name

        return d

    def get_metadata(self, pattern, fields, concatenate=True):
        '''Get metadata'''
        fmt = pattern.split('.')[-1]

        col_dict = self.get_columns_format(
                fmt,
                fields=fields)
        exp_fns = fnmatch.filter(self._zf.namelist(), pattern)

        metadata = []
        for fn in exp_fns:
            with self._zf.open(fn, 'r') as f:
                table = pd.read_csv(
                        io.StringIO(f.read().decode()),
                        sep='\t',
                        header=None,
                        usecols=col_dict.keys())
                table = (table.rename(columns=col_dict)
                              .set_index('Id')
                              .sort_values('CreateDate'))
                metadata.append(table)
        if concatenate:
            metadata = pd.concat(metadata)
        return metadata

    @staticmethod
    def convert_data_for_google_sheet(data):
        cols = ['Plate.barcode', 'dNTP.batch', 'oligodT.order.no', 'plate.type',
                'preparation.site', 'date.prepared', 'data.sorted',
                'tissue', 'subtissue', 'Mouse ID (age_#_sex)',
                'FACS.selection', 'nozzle.size', 'FACS.instument',
                'Experiment ID']

        data = data.copy()

        data['dNTP.batch'] = ''
        data['FACS.selection'] = ''
        data['preparation.site'] = 'Biohub'
        data['date.prepared'] = ''
        data['oligodT.order.no'] = ''
        data.rename(columns={
            'PlateBarcode': 'Plate.barcode',
            'PlateType': 'plate.type',
            'CreateDate': 'data.sorted',
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
            data.at[ix, 'data.sorted'] = datum['data.sorted'].split()[0].replace('-', '')[2:]

        return data[cols]


# Script
if __name__ == '__main__':


    mice_pairs = ('1_6_M 1_7_M', '3_8_M 3_9_M', '3_10_M 3_11_M')
    data = []
    for mice in mice_pairs:
        e = SonyExpdat(mice)
        e.get_filename()
        datum = e.get_full_metadata()
        data.append(datum)
    data = pd.concat(data)
    data_gs = e.convert_data_for_google_sheet(data)

    data_gs.to_csv('../data/tmp/week2clark.tsv', sep='\t', header=False, index=False)
