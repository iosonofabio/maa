# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/06/17
content:    Classes and functions to support SONY's SH800 expdat format.
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

    def __repr__(self):
        return self.__class__.__name__+'(mice='+self.mice+', cytometer='+self.cytometer+')'

    def get_filename(self):
        import glob
        fdn_expdat = '../data/cytometer files/'
        return glob.glob(fdn_expdat+'*'+self.mice+'/'+'*'+self.mice+'*'+self.cytometer+'.expdat')[0]

    def get_metadata(self):
        with zipfile.ZipFile(self.get_filename(), 'r') as zf:
            self._zf = zf

            # Let's get the experiments first (top down)
            exps = self.get_raw_metadata(
                    'Experiments/*.ex',
                    fields=('Id', 'Name', 'Investigator', 'CreateDate'))

            # Then we get the sample groups
            sample_groups = self.get_raw_metadata(
                    'SampleGroups/*.sg',
                    fields=('Id', 'Name', 'CreateDate', 'Experiment_Id'))

            # Then we get the tubes
            tubes = self.get_raw_metadata(
                    'Tubes/*.t',
                    fields=('Id', 'Name', 'CreateDate', 'SampleGroup_Id'))


            # Then, we get the data sources that could have the plate barcode
            data_sources = self.get_raw_metadata(
                    'DataSources/*.ds',
                    fields=('Id', 'Name', 'CreateDate', 'SortResult_Id', 'SortSetting_Id', 'Tube_Id'))

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

            # Find plate barcode
            platename = data_source['Name'].replace(' ', '').rstrip('\\')
            platetype = ''

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
            exname_mice = re.findall('\d+_\d+_[MF](?:_P\d)?', exname)
            tubename_mice = re.findall('\d+_\d+_[MF](?:_P\d)?', data_source['Tube_Name'])
            if len(tubename_mice) == 1:
                mouse_name = tubename_mice[0]
            elif len(exname_mice) == 1:
                mouse_name = exname_mice[0]
            elif ('Test' in data_source['Tube_Name']) or ('trash' in data_source['Tube_Name']):
                mouse_name = 'None_0'
            else:
                raise ValueError('Mouse name not found for this plate: '+str(data_source))

            # The ExpId has to deal with the 2 mice/day policy shift
            mouse_n = int(mouse_name.split('_')[1])
            if mouse_n <= 5:
                exp_id = 'exp'+str(1 + mouse_n)
            else:
                exp_id = 'exp'+str(7 + (mouse_n - 6) // 2)
                # Mice after 5 are parabiosis, but I sometimes forgot the pair
                if '_P' not in mouse_name:
                    mouse_name += '_P'+str(1 + (mouse_n - 12) // 2)

            data_sources.at[ix, 'PlateBarcode'] = platename
            data_sources.at[ix, 'Mouse'] = mouse_name
            data_sources.at[ix, 'Tissue'] = tissue
            data_sources.at[ix, 'Subtissue'] = subtissue
            data_sources.at[ix, 'PlateType'] = platetype
            data_sources.at[ix, 'ExpId'] = exp_id
            #data_sources.at[ix, 'BarcodeFull'] = platename+' '+mouse_name+' '+tissue+' '+subtissue

        # These are fixed at Clark
        if self.cytometer == 'Clark':
            data_sources['Nozzle'] = '100'
            data_sources['CytometerLocation'] = self.cytometer
        else:
            raise ValueError('Cytometer is not Clark, nozzle size??')

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

    def get_raw_metadata(self, pattern, fields, concatenate=True):
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
                table = table.rename(columns=col_dict).set_index('Id')

                if 'CreateDate' in table.columns:
                    table.sort_values('CreateDate', inplace=True)
                metadata.append(table)
        if concatenate:
            metadata = pd.concat(metadata)
        return metadata

    def get_index_metadata(self):
        with zipfile.ZipFile(self.get_filename(), 'r') as zf:
            self._zf = zf

            # Let's get the experiments first (top down)
            sort_settings = self.get_raw_metadata(
                    'SortSettings/*.ss',
                    fields=('Id', 'Index', 'SortType', 'Tube_Id'))

            # Then we get the sample groups
            well_sort_settings = self.get_raw_metadata(
                    'WellSortSettings/*.wss',
                    fields=('Id', 'Row', 'Column', 'SortGate', 'SortSetting_Id'))

            del self._zf


        well_sort_settings['Well'] = well_sort_settings.apply(lambda x: chr(65 + x['Row'])+str(x['Column']+1), axis=1)

        # Keep only wells that have a clean record of metadata
        meta = self.get_metadata()
        meta['DataSource_Id'] = meta.index
        meta['SortSetting_Index'] = sort_settings.loc[meta['SortSetting_Id'].values, 'Index'].values
        meta['SortType'] = sort_settings.loc[meta['SortSetting_Id'].values, 'SortType'].values

        well_sort_settings = well_sort_settings.loc[well_sort_settings['SortSetting_Id'].isin(meta['SortSetting_Id'])]

        meta.set_index('SortSetting_Id', inplace=True, drop=False)
        for col in ('Mouse', 'Tissue', 'Subtissue', 'PlateBarcode',
                    'ExpId', 'Nozzle', 'CytometerLocation', 'Tube_Name',
                    'SortType'):
            well_sort_settings[col] = well_sort_settings.apply(
                    lambda x: meta.at[x['SortSetting_Id'], col], axis=1)

        return well_sort_settings
