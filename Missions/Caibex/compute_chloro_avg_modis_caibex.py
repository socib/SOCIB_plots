#!/usr/bin/env python
#
# compute_chloro_avg_modis_caibex.py
#
#
#------------------------------------------------------------------------------

import os
import glob 
import numpy as np
import math
import netCDF4 as netcdf
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy import stats                     # to compute "nanmean", abstent in my numpy

# Clean 
os.system('clear')

dowrite=1

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'h'

# File and directory names
datadir='/data_local/Satellite/MODIS/data/daily/chlorophyll/'
databasename='T2009*.nc'
outputdir='/data_local/Satellite/MODIS/data/daily/avg/'
outputfile='chloro_avg_Terra_Caibex.nc'

valex=-999.0
#------------------------------------------------------------------------------------

# Loop on files
filelist = sorted(glob.glob(datadir+databasename))
nfiles = len(filelist)

# Load coordinates in the first file
file2load=filelist[0]
nc=netcdf.Dataset(file2load)
lon0=nc.variables['lon'][:]
lat0=nc.variables['lat'][:]
ttime=nc.variables['time'][:]
nc.close()

# Allocate
chloromatrix=np.zeros((nfiles,len(lat0),len(lon0)))

tindex=0
for chlorofile in filelist:

    # Load data from current file
    nc=netcdf.Dataset(chlorofile)
    field = nc.variables['chloro'][:]
    ttime = ttime=nc.variables['time'][:]
    nc.close()
    field = field.squeeze()
    #field=np.flipud(field)

    chloromatrix[tindex]=field
    tindex+=1

chloromatrix2=np.ma.masked_where(np.isnan(chloromatrix),chloromatrix)
chloroavg=np.mean(chloromatrix2,axis=0)

if dowrite==1:
    print 'Write into new NetCDF file'
    
    rootgrp = Dataset(outputdir+outputfile,'w', format='NETCDF4')
    #print rootgrp.file_format

    rootgrp.author = 'ctroupin'
    rootgrp.institute = 'SOCIB'
     
    # Create dimensions
    lat = rootgrp.createDimension('lat', len(lat0))
    lon = rootgrp.createDimension('lon', len(lon0))

    # Create variables
    latitudes = rootgrp.createVariable('latitude','f4',('lat'))
    longitudes = rootgrp.createVariable('longitude','f4',('lon'))
    chloro = rootgrp.createVariable('chloro','f4',('lat','lon',))

    # Put attributes
    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    chloro.units = 'degrees celsius'
    chloro.long_name = 'Sea Surface Temperature'
    chloro.missing_value = valex

    # write variables
    latitudes[:]= lat0
    longitudes[:] = lon0
    chloro[:,:] = chloroavg
    
    # close NetCDF
    rootgrp.close()

    

