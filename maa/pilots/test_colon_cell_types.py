# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/11/17
content:    Try to see where in the sorting plots are successful and failed
            cells for different colon cell types (after RNA-Seq annotation).
'''
# Modules
import os
import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


# Functions
def parse_facs_plate(plate):
    import glob
    import fcsparser

    # FIXME @Joan: change this!!
    fdn = '../../data/index_fcs_colon/'
    fn_fcs = glob.glob(fdn+'*'+plate+'*.fcs')[0]
    fn_index = glob.glob(fdn+'*'+plate+'*_Index.csv')[0]

    out = {}

    meta, data = fcsparser.parse(fn_fcs, reformat_meta=True)
    out['fcs_meta'] = meta
    out['fcs_data'] = data

    data_index = pd.read_csv(fn_index, sep=',', index_col='Index')
    data_index.index = pd.Index(
            [plate+'_'+i for i in data_index.index],
            name='name')
    out['index_data'] = data_index

    # Post-processing:
    # - index data may use BSC instead of SSC
    if 'BSC-A' in data_index.columns:
        data_index.rename(
                columns={
                    'BSC-A': 'SSC-A',
                    'BSC-H': 'SSC-H',
                    'BSC-W': 'SSC-W'},
                inplace=True)

    # - index data may use FITC-Compensated, rename to FITC
    rename_comp = {}
    for x in data_index.columns:
        if x.endswith('-Compensated'):
            rename_comp[x] = x[:-len('-Compensated')]
    data_index.rename(columns=rename_comp, inplace=True)

    # - index data is compensated, but FCS is not by default
    if '$COMP' in out['fcs_meta']:
        com = list(map(float, out['fcs_meta']['$COMP'].split(',')))[1:]
        n_chan = int(np.sqrt(len(com)))
        # This is the right order, Fortran style I guess
        com = np.reshape(com, (n_chan, n_chan)).T
        fluors = out['fcs_meta']['_channel_names_'][-n_chan:]
        com = pd.DataFrame(
                data=com,
                index=fluors,
                columns=fluors)
        out['fcs_data'].loc[:, fluors] = out['fcs_data'].loc[:, fluors].values @ com.values

    # - add plate name to all FCS lines
    out['fcs_data'].loc[:, 'plate'] = plate

    return out

def parse_facs_map():
    #fn = '../../data/index_fcs_colon/FACSmap.csv'
    fn = '../../data/index_fcs_colon/FACSmap_Taichi.csv'
    out = pd.read_csv(fn, sep=',', index_col=0)

    # FIXME: use actual annotations instead of Bob's clusters
    out.rename(columns={'cluster': 'cell_type_call'}, inplace=True)
    out.loc[:, 'cell_type_call'] = out.loc[:, 'cell_type_call'].astype('U30')

    # NOTE: plate 937 is probably a misannotation for 585
    out.index = [i if 'MAA000937' not in i else i.split('.')[0]+'.'+'MAA000585'
            for i in out.index]
    out['plate'] = [i.split('.')[1] for i in out.index]

    out.index.name = 'id'
    out['name'] = ['{1}_{0}'.format(*i.split('.')) for i in out.index]

    out.set_index('name', drop=True, inplace=True)
    return out


# Script
if __name__ == '__main__':

    dead_stain = 'CD45-DAPI: Pacific Blue-A'

    facs_map = parse_facs_map()
    plates = np.unique(facs_map['plate'])
    subtissues = {}
    for plate, datum in facs_map.groupby('plate'):
        subtissues[plate] = np.unique(datum['subtissue'])[0]

    facs_data = {}
    for ip, plate in enumerate(plates):
        print('n. {:}, {:}'.format(ip+1, plate))

        facs_datum = parse_facs_plate(plate)
        ind_bool = facs_datum['index_data'].index.isin(facs_map.index)
        ind = facs_datum['index_data'].index[ind_bool]
        facs_datum['index_data']['Good quality'] = ind_bool

        facs_datum['index_data']["cell type call"] = "missing"
        facs_datum['index_data'].loc[ind, "cell type call"] = \
            facs_map.loc[ind, "cell_type_call"]

        facs_datum['index_data']["n reads"] = 0
        facs_datum['index_data'].loc[ind, "n reads"] = facs_map.loc[ind, "Reads"]

        facs_datum['index_data']["n genes"] = 0
        facs_datum['index_data'].loc[ind, "n genes"] = facs_map.loc[ind, "Genes"]

        facs_data[plate] = facs_datum

        #for ig, groupby in enumerate(('Good quality', 'cell type call')):
        for ig, groupby in enumerate(('cell type call',)):
            n_groups = len(np.unique(facs_datum['index_data'][groupby]))

            fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 9))
            axs = axs.ravel()
            colors = sns.color_palette('husl', n_colors=n_groups)
            ax = axs[0]
            for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                datum.plot(
                        kind='scatter',
                        x='FSC-A',
                        y='SSC-A',
                        label=qual,
                        ax=ax,
                        color=colors[i],
                        zorder=5 - 0.5 * (qual == 'missing'),
                        alpha=0.7 - 0.65 * (qual == 'missing'),
                        )

            facs_datum['fcs_data'].plot(
                    kind='scatter',
                    x='FSC-A',
                    y='SSC-A',
                    ax=ax,
                    logy=True,
                    xlim=(1e2, 1e7),
                    ylim=(1e3, 1e7),
                    color=[0.4] * 3,
                    grid=True,
                    zorder=4,
                    alpha=0.4,
                    )

            ax.legend(loc='lower right', title=groupby)

            ax = axs[1]
            for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                datum.plot(
                        kind='scatter',
                        x='FSC-A',
                        y='SSC-A',
                        ax=ax,
                        color=colors[i],
                        zorder=5 - 0.5 * (qual == 'missing'),
                        alpha=0.7 - 0.65 * (qual == 'missing'),
                        )

            facs_datum['fcs_data'].plot(
                    kind='scatter',
                    x='FSC-A',
                    y='SSC-A',
                    ax=ax,
                    logy=True,
                    xlim=(1e2, 1e6),
                    ylim=(1e3, 1e6),
                    color=[0.4] * 3,
                    grid=True,
                    zorder=4,
                    alpha=0.4,
                    )

            ax = axs[2]
            for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                datum.plot(
                        kind='scatter',
                        x='TIME',
                        y=dead_stain,
                        ax=ax,
                        color=colors[i],
                        zorder=5 - 0.5 * (qual == 'missing'),
                        alpha=0.7 - 0.65 * (qual == 'missing'),
                        )

            facs_datum['fcs_data'].plot(
                    kind='scatter',
                    x='TIME',
                    y=dead_stain,
                    ax=ax,
                    logy=True,
                    ylim=(1e1, 1e6),
                    color=[0.4] * 3,
                    grid=True,
                    zorder=4,
                    alpha=0.4,
                    )

            #ax = axs[3]
            #for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
            #    datum.plot(
            #            kind='scatter',
            #            x='FSC-A',
            #            y=dead_stain,
            #            ax=ax,
            #            color=colors[i],
            #            zorder=5 - 0.5 * (qual == 'missing'),
            #            alpha=0.7,
            #            )

            #facs_datum['fcs_data'].plot(
            #        kind='scatter',
            #        x='FSC-A',
            #        y=dead_stain,
            #        ax=ax,
            #        logy=True,
            #        xlim=(1e2, 1e6),
            #        ylim=(1e1, 1e6),
            #        color=[0.4] * 3,
            #        grid=True,
            #        zorder=4,
            #        alpha=0.4,
            #        )

            posx = 'CD44: APC-A'
            posy = 'CD66a: PE-A'

            ax = axs[3]
            for i, (qual, datum) in enumerate(facs_datum['index_data'].groupby(groupby)):
                datum.plot(
                        kind='scatter',
                        x=posx,
                        y=posy,
                        ax=ax,
                        color=colors[i],
                        zorder=5 - 0.5 * (qual == 'missing'),
                        alpha=0.7 - 0.65 * (qual == 'missing'),
                        )

            facs_datum['fcs_data'].plot(
                    kind='scatter',
                    x=posx,
                    y=posy,
                    ax=ax,
                    logx=True,
                    logy=True,
                    xlim=(1e2, 1e5),
                    ylim=(1e1, 1e6),
                    color='grey',
                    edgecolor='none',
                    grid=True,
                    zorder=4,
                    alpha=0.03,
                    s=10,
                    )

            fig.suptitle(plate+', '+subtissues[plate])
            plt.tight_layout(rect=(0, 0, 1, 0.96))

            fig.savefig('../../figures/{:}_{:}.png'.format(plate, str(ig+1)))
            plt.close(fig)

        ##FIXME
        #if ip == 1:
        #    break

    #plt.ion()
    #plt.show()
