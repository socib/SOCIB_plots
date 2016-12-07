#!/usr/bin/env python

# diva_plot_results.py
#
# Plot the gridded field corresponding to the analysis or the error
# (netCDF file)
#
# For error plotting: change line
# field = nc.variables['analyzed_field'][:]
# into
# field = nc.variables['error_field'][:]
#
# http://modb.oce.ulg.ac.be/mediawiki/index.php/Diva_python
#------------------------------------------------------------------------------

import os, glob
import numpy as np
import math
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors

varname='temp'
datadir='/home/ctroupin/DataOceano/CAIBEX_campaign/data/Thermosal/tss_all/'
datafile='thermosal_sst2.dat'
resultdir='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/results/Horizontal/V2/'
resultfile='thermosal.nc'
figdir='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/figures/results/horizontal/V2/'

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'l'
maxerror=0.25    # Maximal relative error for the plot


# region of interest
lonmin = -12.25 
lonmax = -9.5 
latmin = 29.75
latmax = 31.75

# Limits in the matrices (determined emperically)
imin,imax,jmin,jmax=450,800,1500,1800

cmap = plt.cm.spectral

# Compute min and max values 
vmin=17.0
vmax=23.5
norm = colors.Normalize(vmin=vmin,vmax=vmax)
nlevels=150
levels2plot=np.arange(vmin,vmax+0.0001,0.5)
newticks=np.arange(np.floor(vmin),np.ceil(vmax)+0.0001,1.0)
newticks=newticks.round(2)


levels2plot=np.arange(12,23,0.5)

# Where to put the lon/lat ticks
dlon,dlat=.5,.5
lon2plot = np.arange(-12.5,-9.,dlon)
lat2plot = np.arange(29.5,32.,dlon)

# Spacing between x/y labels
dlon = 0.5
dlat = 0.5

# prepare projection
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                    urcrnrlon=lonmax,urcrnrlat=latmax,  \
                    lat_ts=0.5*(latmin+latmax),\
                    resolution=basemap_resolution)

#------------------------------------------------------------------------------------

# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)


# Loop on result files

# Figure name
figbasename = os.path.basename(resultfile).split('.')[0]
print figbasename
num=figbasename[-5:]


londata1,latdata1=np.loadtxt(datadir+datafile,usecols=(0,1),unpack=True)

# Load results from file
nc=netcdf.Dataset(resultdir+resultfile)
lon = nc.variables['x'][:]
lat = nc.variables['y'][:]
field = nc.variables['analyzed_field'][:]
error = nc.variables['error_field'][:]
nc.close()

# Mask land values
valex=field.min()
field = np.ma.masked_array(field,field==valex)
##    field = np.ma.masked_where(error>=maxerror,field)

# Make the plot
fig=plt.figure()
ax = fig.add_subplot(111)

m.ax=ax

llon,llat=np.meshgrid(lon,lat)
                             
x,y = m(llon,llat)
x1,y1 = m(londata1,latdata1)

contour=m.contour(x,y,field,levels2plot,cmap=cmap,zorder=2)
plt.clabel(contour,levels2plot[0::2],inline=1,fmt='%1.1f',fontsize=14,fontname='Times New Roman')

#plt.plot(x1,y1,'k--')

# Add grid, coastline and continent
m.drawcoastlines(ax=ax,zorder=4)
m.fillcontinents(color='0.7', ax=ax,zorder=3)
m.drawmeridians(lon2plot,labels=[0, 0, 0, 1],fontsize=18,fontname='Times New Roman')
m.drawparallels(lat2plot,labels=[1, 0, 0, 0],fontsize=18,fontname='Times New Roman')

# Export figure and display it
##    plt.savefig(figdir+figbasename, dpi=300, facecolor='w', edgecolor='w',
##                 transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
plt.close()

