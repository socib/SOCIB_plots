#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import netCDF4 as netcdf
import numpy as np
import cf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

datadir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/"
datafilelist = sorted(glob.glob(datadir + "*.nc"))

coordinates = (-6.25, 42., 30., 48.)

nfiles = len(datafilelist)
print str(nfiles) + ' files'

for datafiles in datafilelist:
    print datafiles

    #f = cf.read(datafiles)
    with netcdf.Dataset(datafiles, 'r') as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]

    plt.plot(lon, lat, 'ko', ms=1)

plt.show()
plt.close()