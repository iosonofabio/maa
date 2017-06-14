# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/06/17
content:    Filenames factory.
'''
# Modules
import os
import glob


# Globals
if 'MAA_ROOT_DATA_FOLDER' in os.environ:
    maa_root_data_folder = og.getenv('MAA_ROOT_DATA_FOLDER').rstrip('/')+'/'
else:
    maa_root_data_folder = '/home/fabio/university/postdoc/mouse ageing atlas/maa/data/'

cytometry_data_folder = maa_root_data_folder+'cytometer files/'


# Functions
def get_fcs_filenames(mouse, tissue):
    '''Get FCS filenames from one sort'''
    fdn = cytometry_data_folder+mouse+'/'+tissue+'/'
    # FIXME: we should SORT the files!
    fns = glob.glob(fdn+'*.fcs')
    # Muscle and other tissues have FCS files for index sorts as well
    fns = [fn for fn in fns if '_INX_' not in fn]
    return fns


def get_fcs_filename(plate):
    '''Get FCS filename from one specific plate'''
    fdn = cytometry_data_folder
    fns = glob.glob(fdn+'**/*'+plate+'*.fcs', recursive=True)
    # Muscle and other tissues have FCS files for index sorts as well
    fns = [fn for fn in fns if '_INX_' not in fn]
    if len(fns) == 0:
        raise IOError('File not found: plate '+plate)
    elif len(fns) > 1:
        raise IOError('More than one file found: plate '+plate+
                      ' found: '+', '.join(fns))
    return fns[0]


def get_index_sort_filenames(mouse, tissue):
    '''Get FCS filenames from one sort'''
    fdn = cytometry_data_folder+mouse+'/'+tissue+'/'
    # Check other tissues on Arias
    if tissue in ('muscle',):
        fns = glob.glob(fdn+mouse+'_INX_*.fcs')
    else:
        # FIXME: we should SORT the files!
        fns = glob.glob(fdn+'*_Index.csv')
    return fns


def get_index_sort_filename(plate):
    '''Get FCS filename from one specific plate

    NOTE: sometimes sorts get adjusted in the middle of a plate, so we CAN have
    more than one file per plate here
    '''
    fdn = cytometry_data_folder
    # Muscle and other tissues have FCS files for index sorts as well
    fns = glob.glob(fdn+'**/*'+plate+'*_Index.csv', recursive=True)
    if len(fns) == 0:
        fns = glob.glob(fdn+'**/*_INX_*'+plate+'.fcs', recursive=True)
    if len(fns) == 0:
        raise IOError('File not found: plate '+plate)
    elif len(fns) > 1:
        return fns
    return fns[0]
