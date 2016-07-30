#!/usr/bin/env python
#
# -------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
import datetime, time, calendar
import scipy.io
import matplotlib.text as text
from osgeo import gdal

from alborex_functions import *

dlon, dlat = 1.0, 1.0
coordinates = np.array((-5, 2., 35., 39.))
res = 'i'

stationfile = '/home/ctroupin/SOCIB/Facilities/Vessel/20151126_ALMO/stations_list.txt'
datadir = '/home/ctroupin/DataOceano/Satellite/Chlorophyll/ALMO201511/'
figdir = '/home/ctroupin/DataOceano/Satellite/Chlorophyll/ALMO201511/figures/'
chlorofilelist = sorted(glob.glob(datadir + '*nc'))
nfiles = len(chlorofilelist)
valex = 999

# Read point list from file
latpts, lonpts = np.loadtxt(stationfile, unpack='True', usecols=(0,1))

# Colormap
cmapchloro = plt.cm.YlGnBu_r

# Compute min and max values 
chloromin, chloromax = 0.05, 1.
normchloro = colors.Normalize(vmin=chloromin, vmax=chloromax)
boundchloro = np.arange(chloromin, chloromax + .001, 1)


fig, m, ax = prepare_map(coordinates, res)
plt.close()

lonpts, latpts = m(lonpts, latpts)

for f in range(0, nfiles):
    # Load L2 data
    fig, m, ax = prepare_map(coordinates, res)
    figname = os.path.basename(chlorofilelist[f])[:-3].replace('.', '_')

    with netCDF4.Dataset(chlorofilelist[f], 'r') as nc:
        lon = nc.groups['navigation_data'].variables['longitude'][:]
        lat = nc.groups['navigation_data'].variables['latitude'][:]
        chla = nc.groups['geophysical_data'].variables['chlor_a'][:]
        timechla = nc.time_coverage_start

        NN = 1
        # np.ma.masked_inside(chla, 0, 10)

        # x,y=m(lon[latstart:-latend:NN,lonstart:-lonend:NN],lat[latstart:-latend:NN,lonstart:-lonend:NN])
        x, y = m(lon, lat)
        chloropcm = m.pcolormesh(x, y, chla, cmap=cmapchloro, norm=normchloro, edgecolor='none')

        cbar = fig.colorbar(chloropcm, cmap=cmapchloro, norm=normchloro, orientation='vertical', pad=0.025, aspect=15,
                            shrink=1, extend='both')

        m.plot(lonpts, latpts, 'ko', ms=2, markerfacecolor='k', alpha=0.75)

        m.drawparallels(np.arange(coordinates[2], coordinates[3], dlat), linewidth=0.5,
                        labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
        m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.5,
                        labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

        m.fillcontinents(ax=ax, color='0.5', zorder=2)
        m.drawcoastlines(ax=ax, linewidth=0.2)


        plt.title(timechla, fontsize=20)

        plt.savefig(figdir + figname + '_b', dpi = 300, facecolor = 'w', edgecolor = 'w',
                    transparent = False, bbox_inches = 'tight', pad_inches = 0.1)

#        plt.show()
        plt.close()
