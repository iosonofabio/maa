# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/06/17
content:    FACS sorting related plumbing.
'''
# Modules
import os
import numpy as np
import pandas as pd
import xarray as xr


# Classes / functions
class FacsSample:
    def __init__(self, mouse, tissue):
        self.mouse = mouse
        self.tissue = tissue
        self.sorts = [FacsSort(pn) for pn in self.get_platenames()]

    def __repr__(self):
        return 'FacsSample(\''+self.mouse+'\', \''+self.tissue+'\')'

    def __str__(self):
        return 'FACS sample of mouse '+self.mouse+', tissue '+self.tissue

    def get_platenames(self):
        '''Get the names of plates from this sample'''
        from .filenames import get_fcs_filenames
        fns = get_fcs_filenames(self.mouse, self.tissue)
        plates = []
        for fn in fns:
            fnbn = os.path.basename(fn)
            if 'MAA' not in fnbn:
                raise ValueError('No plate found in filename: '+fnbn)
            ind = fnbn.find('MAA')
            plate = fnbn[ind: ind+9]
            plates.append(plate)
        return plates


class FacsSort:
    def __init__(self, plate):
        self.plate = plate

    def __repr__(self):
        return 'FacsSort(\''+self.plate+'\')'

    def __str__(self):
        return 'FACS sort of plate '+self.plate

    def _parse_fcs(self):
        from .filenames import get_fcs_filename
        import fcsparser
        if not hasattr(self, 'fn'):
            self.fn = get_fcs_filename(self.plate)
        self._fcs_metadata, self._fcs_data = fcsparser.parse(
                self.fn,
                meta_data_only=False,
                reformat_meta=True)

    def get_fcs_data(self):
        '''Get FCS data from this plate'''
        if not hasattr(self, '_fcs_data'):
            self._parse_fcs()
        return self._fcs_data.copy()

    def get_fcs_metadata(self):
        '''Get FCS data from this plate'''
        if not hasattr(self, '_fcs_metadata'):
            self._parse_fcs()
        return self._fcs_metadata.copy()

    def get_channels(self):
        '''Get FACS channels'''
        return self.get_fcs_metadata()['_channels_']

    def plot_fcs_data(
            self,
            axes=(('FSC-A', 'SSC-A'),),
            scales=(('linear', 'linear'),),
            ):
        '''Scatter plot the FCS data'''
        import matplotlib.pyplot as plt
        import seaborn as sns

        if len(axes) != len(scales):
            raise ValueError('axes and scales must have the same length')
        nplots = len(axes)

        if nplots == 1:
            nrows, ncols = 1, 1
        elif nplots == 2:
            nrows, ncols = 1, 2
        elif nplots == 3:
            nrows, ncols = 1, 3
        elif nplots == 4:
            nrows, ncols = 2, 2
        elif nplots in (5, 6):
            nrows, ncols = 2, 3
        elif nplots in (7, 8):
            nrows, ncols = 2, 4
        elif nplots == 9:
            nrows, ncols = 3, 3
        elif nplots in (10, 11, 12):
            nrows, ncols = 3, 4
        elif nplots in (13, 14, 15):
            nrows, ncols = 3, 5
        elif nplots == 16:
            nrows, ncols = 4, 4
        else:
            nrows = int(np.floor(np.sqrt(nplots)))
            if nrows**2 == nplots:
                ncols = nrows
            else:
                ncols = nrows + 1

        width = 4 * ncols
        height = 3.3 * nrows

        fs = 14
        fig, axs = plt.subplots(
                nrows, ncols,
                sharex=False, sharey=False,
                figsize=(width, height))
        if nplots == 1:
            axs = [axs]
        else:
            axs = axs.ravel()

        data = self.get_fcs_data()
        for (ax, (xaxis, yaxis), (xscale, yscale)) in zip(axs, axes, scales):
            ax.scatter(data[xaxis], data[yaxis],
                       s=30, color='b', alpha=0.4)

            if xscale == 'log':
                ax.set_xlim(xmin=1e-1)
                ax.set_xscale(xscale)
            if yscale == 'log':
                ax.set_ylim(ymin=1e-1)
                ax.set_yscale(yscale)

            ax.set_xlabel(xaxis, fontsize=fs)
            ax.set_ylabel(yaxis, fontsize=fs)

        fig.suptitle(self.plate, fontsize=fs)

        plt.tight_layout(rect=(0, 0, 1, 0.97))

        return {'fig': fig, 'axs': axs}
