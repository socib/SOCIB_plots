#!/usr/bin/python

import netCDF4
import glob
import os

datadir = "/home/ctroupin/DataOceano/CMEMS/INSITU_BAL_NRT_OBSERVATIONS_013_032/monthly/mooring/201604/"
datafilelist = sorted(glob.glob(os.path.join(datadir, '*.nc')))

for datafiles in datafilelist:
    print(datafiles)

    
