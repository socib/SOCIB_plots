#!/usr/bin/env python
#
#
# Convert SST file into Diva compatible file format
#
#------------------------------------------------------------------------------

import os
import glob 
import numpy as np
import math
import netCDF4 as netcdf
import matplotlib.pyplot as plt
#import gmtColormap
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors

#-------------
# User options
#-------------

sstmin,sstmax=10.,25.0
plotdaily=1

# Resolution for coastline 
basemap_resolution = 'h'

# File and directory names
datadir='/data_local/Satellite/AVHRR/NAR/SSTonly/4Caibex/'
databasename='200908*'
outputdir=datadir+'ascii/'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'

if not(os.path.exists(outputdir)):
    os.makedirs(outputdir)
    
# Limits in the matrices (determined emperically)
imin,imax,jmin,jmax=450,800,1500,1800


# Loop on files
filelist = sorted(glob.glob(datadir+databasename))
nfiles = len(filelist)

# Load coordinates from coordinate file
nc=netcdf.Dataset(coorddir+coordfile)
lon = nc.variables['lon'][imin:imax,jmin:jmax]
lat = nc.variables['lat'][imin:imax,jmin:jmax]
nc.close()

for SSTfile in filelist:
    filename= os.path.basename(SSTfile)
    figname=filename.split(".")[0]
    figtitle=str(figname[0:4])+'-'+str(figname[4:6])+'-'+str(figname[6:8])

    outputfile=filename.split('.')[0]+'.txt'                                    
    print figname

    # Load data from current file
    nc=netcdf.Dataset(SSTfile)
    field = nc.variables['sst'][0,imin:imax,jmin:jmax]
    nc.close()

    # Select sub-region
    field = field.squeeze()

    # Mask land values
    goodmask = np.where(np.logical_and((field<=sstmax),(field>=sstmin)))

    lon2=lon[goodmask]
    lat2=lat[goodmask]
    field2=field[goodmask]

    array2write=zip(lon2,lat2,field2)
    
    np.savetxt(outputdir+outputfile,array2write)

