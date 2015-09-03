#!/usr/bin/env python
#
# plot all the drufter positions for a given period
# --------------------------------------------------------------------------------

import numpy as np
import netCDF4 as netcdf
import time, datetime
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
import pysocibclient
import os
from alborex_functions import *

doplot = 1  # switch for the plot
cmap = plt.cm.RdYlBu_r

res = 'h'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'  # Coastline
coastfile = 'coatline_cartex.dat'
topodir = '/home/ctroupin/DataOceano/Bathymetry/'
topofile = 'medsea_bathymetry_etopo2.nc'
figdir = '/home/ctroupin/Pictures/SOCIB/'
altimetryfile = "http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/madt/L4/2014/06/nrt_med_allsat_madt_uv_20140625_20140701.nc.gz"
altimetryfile2 = ("http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/" +
                 "madt/L4/2014/06/nrt_med_allsat_madt_h_20140625_20140701.nc.gz")

rcParams['contour.negative_linestyle'] = 'solid'

projectname = 'PERSEUS'

levels2plot = np.arange(-0.3, 0.30001, 0.02)


# region of interest
coordinates = np.array((-6, 10., 35., 45.001))
dlon, dlat = 3.0, 2.0

valex = -999.

tinit = (datetime.datetime(2014, 1, 1, 0, 0, 0)-datetime.datetime(1970,1,1)).total_seconds()
tend = (datetime.datetime(2015, 1, 1, 0, 0, 0)-datetime.datetime(1970,1,1)).total_seconds()

# prepare colormaps
sstmin, sstmax = 18, 29.5
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 2.0)
cmapsst = plt.cm.RdYlBu_r

# Generate lists of platforms
socib_api = pysocibclient.ApiClient()
drifterlist = socib_api.list_platforms(init_datetime="2014-01-00T000000", instrument_type="surface_drifter")

loncoast, latcoast = np.loadtxt(coastdir + coastfile, usecols=(0, 1), unpack=True)
loncoast[loncoast == valex] = np.nan
latcoast[latcoast == valex] = np.nan

# Prepare plot
# ------------

if not os.path.exists(figdir):
    os.makedirs(figdir)

fig, m, ax = prepare_map(coordinates, res)
xc, yc = m(loncoast, latcoast)


ndrifter = len(drifterlist)
print 'Working on ' + str(ndrifter) + ' drifters'
print ' '


ax.set_xlim(coordinates[0], coordinates[1])
ax.set_ylim(coordinates[2], coordinates[3])
m.drawparallels(np.arange(36., 45.0001, dlat), linewidth=0.,
                labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
                labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

m.plot(xc, yc, 'k', lw=0.2, zorder=5)

p1 = m.scatter(0.0, 0.0, c=15, s=0.0, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst)

cbar = fig.colorbar(p1, cmap=cmapsst, norm=normsst, orientation='vertical', pad=0.025, aspect=15, shrink=0.9,
                    extend='both')
cbar.set_ticks(boundsst)
cbar.set_label(r'$^{\circ}$C', fontsize=16, rotation=0, ha='left')

landpoints = np.array([])

for drifter in drifterlist:

    drifter_opendap = drifter.product_list[-1].opendap
    print drifter_opendap
    with netcdf.Dataset(drifter_opendap) as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        gtime = nc.variables['time'][:]

        goodcoord = np.where((lon >= coordinates[0]) & (lon <= coordinates[1]) &
                         (lat >= coordinates[2]) & (lat <= coordinates[3]))[0]

        lon, lat, gtime = lon[goodcoord], lat[goodcoord], gtime[goodcoord]

        goodtime = np.where( np.logical_and( (gtime>=tinit), (gtime<tend)))[0]

        lon = lon[goodtime]
        lat = lat[goodtime]
        # lat = np.ma.masked_outside(gtime, tinit, tend, lat)
        # lon = np.ma.masked_outside(gtime, tinit, tend, lon)

        if len(goodtime)> 0:


            lon, lat = m(lon, lat)
            # pp = map(m.is_land, lon, lat)
            # print sum(pp)
            # landpoints = np.append(landpoints, pp)

            try:
                WTEM = nc.variables['WTEM'][goodcoord]
                WTEM = WTEM[goodtime]

                m.scatter(lon, lat, c=WTEM,
                          s=2, edgecolor='none', cmap=cmapsst, norm=normsst, zorder=6, alpha=0.75)
            except KeyError:
                # print 'No variable WTEM'
                m.plot(lon, lat, 'ko', ms=0.5)

# print 'Total number of points on land : ' + str(sum(landpoints))
# print ' '
# print len(landpoints)

plt.savefig(figdir + 'drifters2014', dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)
# plt.show()
plt.close()
