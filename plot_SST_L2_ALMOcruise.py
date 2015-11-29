#!/usr/bin/env python
# plot_SST_L2.py
# -------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
import datetime, time, calendar

from alborex_functions import *

# File and directory names

dlon, dlat = 1.0, 1.0
coordinates = np.array((-3, 1.5, 35., 39.))
res = 'i'

datadir = "/home/ctroupin/DataOceano/Satellite/SST/ALMOcruise201511/"
figdir = '/home/ctroupin/DataOceano/Satellite/SST/ALMOcruise201511/figures/'
vesselfile = 'http://thredds.socib.es/thredds/dodsC/research_vessel/thermosalinometer/socib_rv-scb_tsl001/L1/2015/11/dep0026_socib-rv_scb-tsl001_L1_2015-11-26.nc'

sstfilelist = sorted(glob.glob(datadir + '*nc'))
nfiles = len(sstfilelist)
valex = 999

# Colormap
cmapsst = plt.cm.RdYlBu_r

# Compute min and max values 
sstmin, sstmax = 16.5, 20.5
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 1.0)

# Load the thermosalinograph data
with netCDF4.Dataset(vesselfile, 'r') as nc:
    lonvessel = nc.variables['LON'][:]
    latvessel = nc.variables['LAT'][:]
    tempvessel = nc.variables['WTR_TEM'][:]

fig, m, ax = prepare_map(coordinates, res)
lonvessel, latvessel = m(lonvessel, latvessel)
plt.close()

for f in range(0, nfiles):
    fig, m, ax = prepare_map(coordinates, res)

    figname = os.path.basename(sstfilelist[f])[:-3].replace('.', '_')

    # Load L2 data



    with netCDF4.Dataset(sstfilelist[f], 'r') as nc:
        lon = nc.groups['navigation_data'].variables['longitude'][:]
        lat = nc.groups['navigation_data'].variables['latitude'][:]
        sst = nc.groups['geophysical_data'].variables['sst'][:]
        qualsst = nc.groups['geophysical_data'].variables['qual_sst'][:]
        ssttime = nc.time_coverage_start

    sst = np.ma.masked_where(qualsst > 0, sst)
    print ssttime

    # Mask
    NN = 1

    x, y = m(lon, lat)
    sstpcm = m.pcolormesh(x, y, sst, cmap=cmapsst, norm=normsst)
    sstscat = m.scatter(lonvessel, latvessel, 15, c=tempvessel,
                        cmap=cmapsst, norm=normsst, edgecolor='None')
    plt.plot(lonvessel, latvessel, 'ko', markersize=1)

    cbar = fig.colorbar(sstpcm, cmap=cmapsst, norm=normsst, orientation='vertical', pad=0.025,
                        aspect=15, shrink=1, extend='both')
    cbar.set_ticks(boundsst)

    m.drawparallels(np.arange(coordinates[2], coordinates[3], dlat), linewidth=0.5,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.5,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

    m.fillcontinents(ax=ax, color='0.5', zorder=2)
    m.drawcoastlines(ax=ax, linewidth=0.2)

    plt.title(ssttime[:10], fontsize=20)

    plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)

#    plt.show()
    plt.close()
