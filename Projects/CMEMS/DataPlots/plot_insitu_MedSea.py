#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import netCDF4 as netcdf
import numpy as np
import cf
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#datadir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/"
datadir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/profiler-glider/201505/"
figdir = "/home/ctroupin/Projects/201501_InsTAC/GliderData/figures/"

datafilelist = sorted(glob.glob(datadir + "*.nc"))

cmapsst = plt.cm.RdYlBu_r
sstmin, sstmax = 16., 21.
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)

coordinates = (-6.25, 30., 30., 45.)
coordinates2 = (-1., 2.5, 38.25, 40.25)

nfiles = len(datafilelist)
print str(nfiles) + ' files'

# m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
#             urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
#             lat_ts=0.5*(coordinates[2]+coordinates[3]),
#             resolution='i')

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
            urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
            lat_ts=0.5*(coordinates[2]+coordinates[3]),
            resolution='i')

fig = plt.figure()
#ax = fig.add_axes([0.1, 0.2, 0.7, 0.5])
ax = fig.add_subplot(111)
m.ax = ax

m.drawcoastlines(ax=ax, linewidth=0.2)
m.fillcontinents(color='0.85', ax=ax)

axins = zoomed_inset_axes(ax, 3, loc=1)

for datafiles in datafilelist[0:10]:
    print datafiles

#   Load the coordinates (not feasible with cf, files are not compliant)
    with netcdf.Dataset(datafiles, 'r') as nc:
        lon0 = nc.variables['LONGITUDE'][:]
        lat0 = nc.variables['LATITUDE'][:]

        lon, lat = m(lon0, lat0)

        print nc.institution
        if 'SOCIB' in nc.institution:
            print '*****************'
            m.plot(lon, lat, 'ko-', ms=1, lw=0.5)
        else:
            m.plot(lon, lat, 'o-', color='0.85', markerfacecolor='0.85', markeredgecolor='0.85', ms=1, lw=0.5)

    m.ax = axins
    m.plot(lon, lat, 'ko-', ms=1, lw=0.5)
    m.ax = ax

# #   Load the temperature
#     f = cf.read(datafiles)
#     try:
#         temperature = f.select("sea_water_temperature")[0].array[0]
#         print " "
#         print len(temperature)
#
#         goodtemp = np.where(temperature > 0)
#         lon, lat, temperature = lon[goodtemp], lat[goodtemp], temperature[goodtemp]
#
#         m.scatter(lon, lat, s=10, c=temperature, edgecolor='none', norm=normsst, cmap=cmapsst)
#     except KeyError:
#         print "No temperature variable"

m.drawcoastlines(ax=axins, linewidth=0.2)
mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

n1, n2 = m(coordinates2[0], coordinates2[2])
n3, n4 = m(coordinates2[1], coordinates2[3])
axins.set_xlim(n1, n3)
axins.set_ylim(n2, n4)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
#plt.savefig(figdir + 'glider201505')
plt.show()
plt.close()
