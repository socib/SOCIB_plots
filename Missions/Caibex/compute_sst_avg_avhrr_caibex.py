#!/usr/bin/env python
#
# compute_sst_avg_avhrr_caibex.py
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
datadir='/data_local/Satellite/AVHRR/NAR/SSTonly/files4avg/'
databasename='2009*.nc'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'
outputdir='/data_local/Satellite/AVHRR/NAR/SSTonly/4Caibex/'
outputfile='SST_avg_Caibex.nc'

# Limits in the matrices (determined emperically)
imin,imax,jmin,jmax=450,800,1500,1800

valex=-999.0
#------------------------------------------------------------------------------------

# Loop on files
filelist = sorted(glob.glob(datadir+databasename))
nfiles = len(filelist)

# Load coordinates from coordinate file
nc=netcdf.Dataset(coorddir+coordfile)
lon0 = nc.variables['lon'][imin:imax,jmin:jmax]
lat0 = nc.variables['lat'][imin:imax,jmin:jmax]
nc.close()

# Allocate
SSTmatrix=np.zeros((nfiles,lon0.shape[0],lon0.shape[1]))

# Loop on the files to compute the average
tindex=0                   
for SSTfile in filelist:
    
    nc=netcdf.Dataset(SSTfile)
    field = nc.variables['sst'][0,imin:imax,jmin:jmax]
    nc.close()
    field[field<0.0]=np.nan
    SSTmatrix[tindex]=field
    tindex+=1
                   
SSTmatrix2=np.ma.masked_where(np.isnan(SSTmatrix ),SSTmatrix)
SSTavg2=np.mean(SSTmatrix2,axis=0)

nlat,nlon=SSTavg2.shape

if dowrite==1:
    print 'Write into new NetCDF file'
    
    rootgrp = Dataset(outputdir+outputfile,'w', format='NETCDF4')
    #print rootgrp.file_format

    rootgrp.author = 'ctroupin'
    rootgrp.institute = 'SOCIB'
     
    # Create dimensions
    lat = rootgrp.createDimension('lat', nlat)
    lon = rootgrp.createDimension('lon', nlon)

    # Create variables
    latitudes = rootgrp.createVariable('latitude','f4',('lat','lon'))
    longitudes = rootgrp.createVariable('longitude','f4',('lat','lon'))
    SST = rootgrp.createVariable('SST','f4',('lat','lon',))

    # Put attributes
    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    SST.units = 'degrees celsius'
    SST.long_name = 'Sea Surface Temperature'
    SST.missing_value = valex

    # write variables
    latitudes[:,:]= lat0
    longitudes[:,:] = lon0
    SST[:,:] = SSTavg2
    
    # close NetCDF
    rootgrp.close()

    

