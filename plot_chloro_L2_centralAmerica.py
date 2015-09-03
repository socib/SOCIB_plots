#!/usr/bin/env python
#
# -------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import datetime, time, calendar
import scipy.io
import matplotlib.text as text
from osgeo import gdal

from alborex_functions import *

# File and directory names

dlon,dlat =  2.0,2.0
coordinates = np.array((-104,-95.,9.5,20))
res = 'l'


chlorodir='/home/ctroupin/'
chlorofile = 'A2002365202500.L2_LAC_OC.nc'
figdir='/home/ctroupin/'

valex=999

# Colormap
cmapchloro=plt.cm.YlGnBu_r

# Compute min and max values 
chloromin,chloromax = 0.05,0.6
normchloro= colors.Normalize(vmin=chloromin,vmax=chloromax)
boundchloro = np.arange(chloromin,chloromax+.001,1)

# newticks  = np.array((0.06,0.1,0.3,1.0,))
newticks  = np.arange(0.0,0.6,0.1)
# newlabels = np.array((0.06,0.1,0.3,1.0,))


fig,m, ax = prepare_map(coordinates,res)

# Load L2 data

with netcdf.Dataset(chlorodir + chlorofile, 'r') as nc:
  lon = nc.groups['navigation_data'].variables['longitude'][:]
  lat = nc.groups['navigation_data'].variables['latitude'][:]
  chla = nc.groups['geophysical_data'].variables['chlor_a'][:]
  timechla = nc.time_coverage_start

# Build date
datechla = timechla[:4]+'-'+timechla[4:6]+'-'+timechla[6:8]

# Mask
NN=1
np.ma.masked_inside(chla,0,10)

#x,y=m(lon[latstart:-latend:NN,lonstart:-lonend:NN],lat[latstart:-latend:NN,lonstart:-lonend:NN])
x,y=m(lon,lat)
chloropcm=m.pcolormesh(x,y,chla,cmap=cmapchloro,norm=normchloro,edgecolor='none')

cbar=fig.colorbar(chloropcm,cmap=cmapchloro,norm=normchloro,orientation='vertical',pad=0.025,aspect=15,shrink=1,extend='both')
cbar.set_ticks(newticks)
#cbar.set_ticklabels(newlabels)

m.drawparallels(np.arange(coordinates[2],coordinates[3],dlat), linewidth=0.5,
                        labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=1)
m.drawmeridians(np.arange(coordinates[0],coordinates[1],dlon), linewidth=0.5,
                        labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16,zorder=1)
m.drawmapscale(-100,17.5,-100,17.5, 200, barstyle='simple', units='km', fontsize=12,zorder=3)

m.fillcontinents(ax=ax,color='w',zorder=2)
m.drawcoastlines(ax=ax)

plt.title(datechla,fontsize=20)

figname=chlorofile[:-3]
plt.savefig(figdir+'chloro20021231', dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)


# plt.show()
plt.close()
