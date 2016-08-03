#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import netCDF4 as netcdf
import numpy as np
import matplotlib.pyplot as plt
import shutil
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


figdir = "/home/ctroupin/Projects/201501_InsTAC/GliderData/figures/"
basedir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/profiler-glider/"
socibdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/profiler-glider/socib_imedea_files/"

coordinates = (-6.25, 30., 30., 45.)
coordinates2 = (-1., 4.0, 38.25, 40.25)

yearlist = (2014, 2015)
monthlist = range(1, 13)

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
            urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
            lat_ts=0.5*(coordinates[2]+coordinates[3]),
            resolution='i')

# Loop on years
for years in yearlist:
    yyyy = str(years)
#   Loop on months
    for months in monthlist:
        mm = str(months).zfill(2)

        datadir = basedir + yyyy + mm + "/"
        datafilelist = sorted(glob.glob(datadir + "*.nc"))
        nfiles = len(datafilelist)
        print str(nfiles) + ' files'

        fig = plt.figure()
#       ax = fig.add_axes([0.1, 0.2, 0.7, 0.5])
        ax = fig.add_subplot(111)
        m.ax = ax

        m.drawcoastlines(ax=ax, linewidth=0.2)
        m.fillcontinents(color='0.85', ax=ax)
        plt.title(yyyy + '-' + mm, fontsize=24)

        axins = zoomed_inset_axes(ax, 3, loc=1)

        print yyyy + mm
        for datafiles in datafilelist[0:10]:

#           Load the coordinates (not feasible with cf, files are not compliant)
            with netcdf.Dataset(datafiles, 'r') as nc:
                lon0 = nc.variables['LONGITUDE'][:]
                lat0 = nc.variables['LATITUDE'][:]

                lon0 = np.ma.masked_outside(lon0, coordinates[0], coordinates[1])
                lat0 = np.ma.masked_outside(lat0, coordinates[2], coordinates[3])

                lon, lat = m(lon0, lat0)

#               print nc.institution
                if 'SOCIB' in nc.institution:
                    print nc.institution
                    print datafiles
                    print " "
#                   Copy the files to a common directory
                    shutil.copy2(datafiles, socibdir)
                    m.ax = ax
                    m.plot(lon, lat, 'ko-', ms=1, lw=0.5)
                    m.ax = axins
                    m.plot(lon, lat, 'ko-', ms=1, lw=0.5)

                elif 'IMEDEA' in nc.institution:
                    print nc.institution
                    print datafiles
                    print(" ")
#                   Copy the files to a common directory
                    shutil.copy2(datafiles, socibdir)
                    m.ax = ax
                    m.plot(lon, lat, 'bo-', markerfacecolor='b', markeredgecolor='b', ms=1, lw=0.5)
                    m.ax = axins
                    m.plot(lon, lat, 'bo-', markerfacecolor='b', markeredgecolor='b', ms=1, lw=0.5)

                else:
                    m.ax = ax
                    m.plot(lon, lat, 'o-', color='0.85', markerfacecolor='0.85', markeredgecolor='0.85', ms=1, lw=0.5)
                    m.ax = axins
                    m.plot(lon, lat, 'o-', color='0.85', markerfacecolor='0.85', markeredgecolor='0.85', ms=1, lw=0.5)

        m.drawcoastlines(ax=axins, linewidth=0.2)
        m.fillcontinents(ax=axins, color='0.85')

        mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

        n1, n2 = m(coordinates2[0], coordinates2[2])
        n3, n4 = m(coordinates2[1], coordinates2[3])
        axins.set_xlim(n1, n3)
        axins.set_ylim(n2, n4)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        m.ax = ax

#        plt.savefig(figdir + 'Medsea_glider' + yyyy + mm)
#       plt.show()
        plt.close()
