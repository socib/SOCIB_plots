#!/usr/bin/python
__author__ = 'ctroupin'

datadir = '/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/profiler-glider/socib_imedea_files/'

import glob
import netCDF4 as netcdf
import matplotlib.pyplot as plt

for datafiles in sorted(glob.glob(datadir+'*nc')):
    print datafiles
    with netcdf.Dataset(datafiles, 'r') as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]

        print 'min. longitude = ' + str(lon.min())
        print 'max. longitude = ' + str(lon.max())
        print 'min. latitude = ' + str(lat.min())
        print 'max. latitude = ' + str(lat.max())
        print ' '


        plt.plot(lon, lat)


plt.show()
plt.close()