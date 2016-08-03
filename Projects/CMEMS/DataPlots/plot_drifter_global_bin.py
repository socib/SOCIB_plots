#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from matplotlib import colors

doplot, dowrite, doread = 1, 0, 1
dobinning = 1

year = 2015
month = 7

figdir = "/home/ctroupin/DataOceano/MyOcean/figures/drifters/"
basedir = "/home/ctroupin/DataOceano/MyOcean/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/" + str(year) + str(month).zfill(2) + '/'
landfile = "/data_local/Bathymetry/ne_10m_land"
file2read = "/home/ctroupin/DataOceano/MyOcean/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/lon_lat_temperature_" + str(year) + str(month).zfill(2) + '.txt'
file2read2 = "/home/ctroupin/DataOceano/MyOcean/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/temperaturemean_" + str(year) + str(month).zfill(2) + '.txt'


sstmin, sstmax = 5., 30.
cmapsst = plt.cm.RdYlBu_r
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1, 2.5)
level2plotsst = np.arange(10, sstmax+.1, 0.01)

dloninterp, dlatinterp = 5., 5.
lon_interp, lat_interp = np.meshgrid(np.arange(-180, 181., dloninterp), np.arange(-70, 71., dlatinterp))


yearlist = (2014, 2015)
monthlist = range(1, 13)

if doplot:
    fig = plt.figure(figsize=(12, 8))
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.00)

    m = Basemap(projection='robin', lon_0=0, resolution='c')
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


if dowrite:
    # Loop on files
    k = 0
    lon, lat, temperature, temperatureQC = np.array([]), np.array([]), np.array([]), np.array([])
    filelist = sorted(glob.glob(basedir+'*.nc'))
    for datafiles in filelist:

        # print datafiles
        with netCDF4.Dataset(datafiles) as nc:

            depth = nc.variables['DEPH'][:]
            # print depth.shape[1]
            if depth.shape[1] == 1:

                try:
                    temperatureQC = nc.variables['TEMP_QC'][:]
                    goodtemperature = np.where(temperatureQC == 1)[0]

                    if len(goodtemperature)>0:
                        temperature = np.append(temperature, nc.variables['TEMP'][goodtemperature])
                        lon = np.append(lon, nc.variables['LONGITUDE'][goodtemperature])
                        lat = np.append(lat, nc.variables['LATITUDE'][goodtemperature])

                except KeyError:
                    # print 'No variable temperature in this file'
                    # temperature = np.nan*np.ones_like(lat)
                    k += 1

    print temperature.shape
    print lon.shape
    print lat.shape
    np.savetxt(file2read, np.c_[lon, lat, temperature], fmt='%3.5f %3.5f %2.2f')


if doread:
    lon, lat, temperature = np.loadtxt(file2read, usecols=(0, 1, 2), unpack=True)
    # temperature = np.ma.masked_outside(temperature, -2.0, 40.)


temperature_mean = np.array([])
if dobinning:

    for lon0, lat0 in zip(lon_interp.flatten(), lat_interp.flatten()):

        goodlon = np.where(np.logical_and((lon < lon0 + dloninterp*0.5), (lon > lon0 - dloninterp*0.5)))
        goodlat = np.where(np.logical_and((lat < lat0 + dlatinterp*0.5), (lat > lat0 - dlatinterp*0.5)))
        goodcoord = np.intersect1d(goodlon, goodlat)

        #print goodcoord
        temperature2calc = temperature[goodcoord]
        #print temperature2calc

        # print len(temperature2calc)
        if len(temperature2calc) > 0:
            temperature_mean = np.append(temperature_mean, temperature2calc.mean())
            # print lon0, lat0, temperature_mean

            x0, y0 = m(lon0, lat0)
            #scat = m.scatter(x0, y0, s=50, c=temperature_mean, marker='s',
            #                 edgecolor='None', cmap=cmapsst, norm=normsst, zorder=4)

        else:
            temperature_mean = np.append(temperature_mean, np.nan)

    np.savetxt(file2read2, temperature_mean)

    temperature_mean = temperature_mean.reshape((lon_interp.shape))
    temperature_mean = np.ma.masked_where(np.isnan(temperature_mean), temperature_mean)

    llon_interp, llat_interp = m(lon_interp, lat_interp)
    pcm = m.pcolormesh(llon_interp, llat_interp, temperature_mean, cmap=cmapsst, norm=normsst)
    # scat = m.scatter(lon, lat, s=2, c=temperature, edgecolor='None', cmap=cmapsst, norm=normsst, zorder=4)

    cbar = plt.colorbar(pcm, cmap=cmapsst, norm=normsst, extend='both', shrink=0.7)
    cbar.set_label('$^{\circ}$C', rotation=0, ha='left')
    plt.title('Temperature from surface drifters\n' + str(year) + '-' + str(month).zfill(2))
    plt.savefig(figdir + 'drifter_binned')
    plt.show()
    plt.close()





