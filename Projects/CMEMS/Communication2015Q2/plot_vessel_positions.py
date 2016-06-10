__author__ = 'ctroupin'

import glob
import netCDF4
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import sys
sys.path.append("/home/ctroupin/Software/Python/seawater-3.3.2")
import seawater

doplot = 0
depthinterp = 10.
figdir = "/home/ctroupin/DataOceano/MyOcean/figures/vessel/"

datadir = "/data_local/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/vessel/201506/"
filelist = sorted(glob.glob(datadir+'*EXRE0173.nc'))

coordinates = np.array((-6.75, 37., 29, 46.))
dlon, dlat = 5., 3.

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
            urcrnrlat=coordinates[3], lat_ts=0.5*(coordinates[2]+coordinates[3]), resolution='i')

if doplot == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m.ax = ax

tempinterp = []

for datafile in filelist:
    print datafile

    with netCDF4.Dataset(datafile) as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]
        temp = nc.variables['TEMP'][:]
        time = nc.variables['TIME']
        time2plot = netCDF4.num2date(time[:], time.units)
        nprofiles = len(lon)
        print 'Number of profiles = ' + str(nprofiles)
        try:
            depth = nc.variables['DEPH'][:]
        except KeyError:
            print "No variable depth"
            try:
                pressure = nc.variables['PRES'][:]
                print pressure.shape
                print "convert pressure to depth"
                depth = np.zeros_like(pressure)
                for ii in range(0, pressure.shape[0]):
                    depth[ii, :] = seawater.eos80.dpth(pressure[ii, :], lat[ii])
            except KeyError:
                print "No variable pressure"


#   Loop on the profiles

    for ii in range(0, nprofiles):
#           Interpolate temperature at given level
#        print 'Depth min. = ' + str(depth[ii,:].min())
#        print 'Depth max. = ' + str(depth[ii,:].max())

        if (depth[ii, :].min()<= depthinterp) & (depth[ii, :].max()>= depthinterp):
            print "Ok to interpolate the profiles"
            gooddepth = np.where(depth[ii, :]>0.)[0]
            tempinterp0 =  np.interp(depthinterp, depth[ii, gooddepth], temp[ii, gooddepth])
            print tempinterp0

            tempinterp.append(tempinterp0)
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # plt.plot(temp[ii, :], -depth[ii,:], 'k-', lw=1)
            # plt.plot(temp[ii,:].mean(), -depthinterp, 'ro', ms=15)
            # plt.plot(tempinterp0, -depthinterp, 'bo', ms=10)
            # plt.show()
            # plt.close()

    lon, lat = m(lon, lat)
    if doplot == 1:
        plt.plot(lon, lat, 'ko', ms=2)


if doplot == 1:
    m.drawcoastlines(ax=ax, linewidth=0.1)
    m.fillcontinents(color = 'gray')
    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

    plt.title('CTD cast locations', fontsize=24)
    plt.show()
    plt.close()