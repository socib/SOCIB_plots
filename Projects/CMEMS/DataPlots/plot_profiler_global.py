#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import cf
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import shutil
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from matplotlib  import colors

doplot = 1

year = 2003
month = 1


figdir = "/home/ctroupin/DataOceano/MyOcean/figures/profiler-glider/"
landfile = "/data_local/Bathymetry/ne_10m_land"

sstmin, sstmax = 5., 30.
cmapsst = plt.cm.RdYlBu_r
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1,2.5)
level2plotsst = np.arange(10, sstmax+.1,0.01)

yearlist = np.arange(2003, 2004)
monthlist = np.arange(1, 7)

m = Basemap(projection='robin',lon_0=0,resolution='c')

for year in yearlist:
    for month in monthlist:

        print 'Working on ' + str(year) + str(month).zfill(2)
        print ' '
        basedir = "/home/ctroupin/DataOceano/MyOcean/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/profiler-glider/" + str(year) + str(month).zfill(2) + '/'

        fig = plt.figure(figsize=(12, 8))
        plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.00)

        shp_info = m.readshapefile(landfile, 'scalerank', drawbounds=True)
        ax = plt.gca()
        ax.cla()

        paths = []
        for line in shp_info[4]._paths:
            paths.append(Path(line.vertices, codes=line.codes))

        coll = PathCollection(paths, linewidths=0, facecolors='grey', zorder=2)

        m = Basemap(projection='robin',lon_0=0,resolution='c')
        # drawing something seems necessary to 'initiate' the map properly
        m.drawcoastlines(color='white', zorder=0)
        ax = plt.gca()
        ax.add_collection(coll)

        # Loop on files
        k, j = 0, 0
        filelist = sorted(glob.glob(basedir+'*.nc'))
        print len(filelist)

        for datafiles in filelist:

            # print datafiles
            with netCDF4.Dataset(datafiles) as nc:
                lon = nc.variables['LONGITUDE'][:]
                lat = nc.variables['LATITUDE'][:]
                #depth = nc.variables['DEPH'][:]
                #distance = haversine( (lat[0], lon[0]), (lat[-1], lon[-1]))

                #print depth.shape[0]
                #print 'Distance = ' + str(distance)
                lon, lat = m(lon, lat)


                plt.plot(lon.mean(), lat.mean(), 'ko', ms=1.5)
                plt.plot(lon, lat, 'bo', ms=0.5)

            j+=1

        plt.savefig(figdir + 'profilers_' + str(year) + str(month).zfill(2))
        # plt.show()
        plt.close()
