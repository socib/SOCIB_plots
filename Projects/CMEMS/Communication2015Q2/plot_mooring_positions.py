__author__ = 'ctroupin'

import glob
import netCDF4
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

doplot = 1
doplotseries = 0
plottimeseries = 0

figdir = "/home/ctroupin/DataOceano/MyOcean/figures/mooring/NRT/"
figdir2 = "/home/ctroupin/Projects/201501_InsTAC/MaterialQ2/"

#mooringdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_TS_REP_OBSERVATIONS_013_041/history/mooring/"
mooringdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/history/mooring/"
mooringlist = sorted(glob.glob(mooringdir+'*.nc'))

coordinates = np.array((-6.75, 37., 29, 46.))
dlon, dlat = 5., 3.

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
            urcrnrlat=coordinates[3], lat_ts=0.5*(coordinates[2]+coordinates[3]), resolution='h')

if doplot == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m.ax = ax

    for moorings in mooringlist:
        print moorings

        with netCDF4.Dataset(moorings) as nc:
            lon = nc.variables['LONGITUDE'][:]
            lat = nc.variables['LATITUDE'][:]
            time = nc.variables['TIME']
            time2plot = netCDF4.num2date(time[:], time.units)
            print
            print len(time)
            print time2plot[0]
            print time2plot[-1]

            if plottimeseries:
                try:
                    temp = nc.variables['TEMP'][:]

                    if doplotseries:
                        fig1 = plt.figure()
                        plt.plot(time2plot, temp)
                        plt.title(str(lat.mean()) + 'N - ' + str(lon.mean()) +'E', fontsize=24)
                        fig1.autofmt_xdate()
                        figname = os.path.basename(moorings).split('.')[0]
                        plt.savefig(figdir + figname)
            #            plt.show()
                        plt.close()

                except:
                    print "No variable temperature in this file"

            print lon.mean()
            print lat.mean()

            lon, lat = m(lon.mean(), lat.mean())
            if doplot == 1:
                plt.plot(lon, lat, 'ko', ms=5)
                if "IR_TS_MO_61198.nc" in moorings:
                    plt.plot(lon, lat, 'bo', ms=5)

    #    f = cf.read(moorings)
    #    temp = f.select('sea_water_temperature')
    #    print len(temp)

if doplot == 1:
    figname = "mooring_map"
    m.drawcoastlines(ax=ax, linewidth=0.1)
    m.fillcontinents(color = 'gray')
    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

    plt.title('Mooring locations', fontsize=24)
    plt.savefig(figdir + figname)
    plt.savefig(figdir2 + figname + ".tif")
    plt.savefig(figdir2 + figname + ".png")
    plt.show()
    plt.close()