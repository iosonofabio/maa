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
    fns = glob.glob(fdn+'*.fcs')
    # FIXME: we should SORT the files!
    return fns


def get_fcs_filename(plate):
    '''Get FCS filename from one specific plate'''
    fdn = cytometry_data_folder
    fns = glob.glob(fdn+'**/*'+plate+'*.fcs', recursive=True)
    if len(fns) == 0:
        raise IOError('File not found: plate '+plate)
    elif len(fns) > 1:
        raise IOError('More than one file found: plate '+plate+
                      ' found: '+', '.join(fns))
    return fns[0]
