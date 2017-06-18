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
import zipfile
import io
import argparse
import numpy as np
import pandas as pd


# Functions / classes
class SonyExpdat:

    def __init__(self, mice):
        self.mice = mice

    def get_filename(self):
        import glob
        fdn_expdat = '../data/cytometer files/'
        return glob.glob(fdn_expdat+'*'+self.mice+'/'+'*'+self.mice+'*Clark.expdat')[0]

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

        # Create new barcode that merges the info
        data_sources['BarcodeFull'] = ''
        for ix, data_source in data_sources.iterrows():
            exname = data_source['Experiment_Name']
            if exname.startswith('MAA '):
                exname = exname[4:]
            platename = data_source['Name'].replace(' ', '')
            data_sources.at[ix, 'BarcodeFull'] = platename+' '+exname

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


# Script
if __name__ == '__main__':


    mice = '3_10_M 3_11_M'
    e = SonyExpdat(mice)
    e.get_filename()
    data = e.get_full_metadata()

    sys.exit()


    # First, expdat is just a ZIP archive (uncompressed)
    fn_expdat = '../data/cytometer files/170616 3_10_M 3_11_M/170616 3_10_M 3_11_M Sony Clark.expdat'
    with zipfile.ZipFile(fn_expdat, 'r') as zf:

        # Let's get the experiments first (top down)
        exps = get_metadata(
                zf, 'Experiments/*.ex',
                fields=('Id', 'Name', 'Investigator', 'CreateDate'))

        # Then we get the sample groups
        sample_groups = get_metadata(
                zf, 'SampleGroups/*.sg',
                fields=('Id', 'Name', 'CreateDate', 'Experiment_Id'))

        # Then we get the tubes
        tubes = get_metadata(
                zf, 'Tubes/*.t',
                fields=('Id', 'Name', 'CreateDate', 'SampleGroup_Id'))


        # Then, we get the data sources that could have the plate barcode
        data_sources = get_metadata(
                zf, 'DataSources/*.ds',
                fields=('Id', 'Name', 'CreateDate', 'SortResult_Id', 'Tube_Id'))

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


        # Create new barcode that merges the info
        data_sources['BarcodeFull'] = ''
        for ix, data_source in data_sources.iterrows():
            exname = data_source['Experiment_Name']
            if exname.startswith('MAA '):
                exname = exname[4:]
            platename = data_source['Name'].replace(' ', '')
            data_sources.at[ix, 'BarcodeFull'] = platename+' '+exname
