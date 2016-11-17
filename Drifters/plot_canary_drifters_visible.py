#!/usr/bin/env python
#
# plot_canary_drifters_visible.py
#
# Plot Pedro's drifers and 
# Add visible image from https://earthdata.nasa.gov/labs/worldview/
# Coastline from http://www.ngdc.noaa.gov/mgg_coastline/
#
# ctroupin, April 2014
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

# GDAL does not use python exceptions by default
gdal.UseExceptions()

doplotbathy = 0
doplotsst =0

# Directory and file names
# ------------------------

imagedir='/home/ctroupin/DataOceano/Satellite/Visible/'
imagefile='nasa-worldview-2014-05-06.tif'
coastdir='/home/ctroupin/DataOceano/Coastlines/'
coastfile='canary_coast.dat'

drifterfile1='http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp029-ieo_svp001/L1/2014/dep0001_drifter-svp029_ieo-svp001_L1_2014-04-20.nc'
drifterfile2='http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp030-ieo_svp002/L1/2014/dep0001_drifter-svp030_ieo-svp002_L1_2014-04-20.nc'
drifterfile3='http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp032-ieo_svp004/L1/2014/dep0001_drifter-svp032_ieo-svp004_L1_2014-04-20.nc'
drifterfile4='http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp033-ieo_svp005/L1/2014/dep0001_drifter-svp033_ieo-svp005_L1_2014-04-20.nc'
drifterfile5='http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp035-ieo_svp007/L1/2014/dep0001_drifter-svp035_ieo-svp007_L1_2014-04-20.nc'
drifterfile6='http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp034-ieo_svp006/L1/2014/dep0001_drifter-svp034_ieo-svp006_L1_2014-04-20.nc'

sstdir='/home/ctroupin/DataOceano/Satellite/SST_dineof/data/Canary/'
sstfile='sst_20140505.nc'
sstcoord='latlon_grid.nc'

figdir='/home/ctroupin/SOCIB/Facilities/Lagrangian/figures/'
figname='drifters_sst_visible2'
figtype1='.eps'
figtype2='.png'
# ----------------------------------------------------------------------------
# region of interest
lonmin=-18.5
lonmax=-13
latmin=25.5
latmax=29.5
dlon=0.5
dlat=0.5

# Create min and max times for the drifter trajectory
# (in seconds after 1970/1/1)
dateinit=calendar.timegm((2013,8,1,0,0,0))
dateend=calendar.timegm((2013,8,15,0,0,0))



# ----------------------------------------------------------------------------
# Read geotiff image
# downloaded from https://earthdata.nasa.gov/labs/worldview/
gtif = gdal.Open(imagedir+imagefile)
gtif.GetProjectionRef()

#Set the plot axis limits to the proper map coordinates.
arr = gtif.ReadAsArray()
trans = gtif.GetGeoTransform()
extent = (trans[0], trans[0] + gtif.RasterXSize*trans[1],
          trans[3] + gtif.RasterYSize*trans[5], trans[3])
##
### Colormap for SST
##cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps/')
##cmapsst = colors.LinearSegmentedColormap('my_colormap',cdict,256)

### Compute min and max values 
vmin=18.5
vmax=21.5
normsst = colors.Normalize(vmin=vmin,vmax=vmax)
nlevels=150
levels2plotsst=np.arange(vmin,vmax+0.0001,(vmax-vmin)/(nlevels-1))
newticks=np.arange(np.floor(vmin),np.ceil(vmax)+0.0001,0.5)
cmapsst = plt.cm.Spectral_r

# ----------------------------------------------------------------------------
# Load extracted coastline
valexc=-999.
lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
lonc=np.ma.masked_where(lonc==valexc,lonc)
latc=np.ma.masked_where(latc==valexc,latc)

# ----------------------------------------------------------------------------
# Load drifter trajectory
nc=netcdf.Dataset(drifterfile1)
londrifter1 = nc.variables['LON'][:]
latdrifter1 = nc.variables['LAT'][:]
nc.close()
nc=netcdf.Dataset(drifterfile2)
londrifter2 = nc.variables['LON'][:]
latdrifter2 = nc.variables['LAT'][:]
nc.close()
nc=netcdf.Dataset(drifterfile3)
londrifter3 = nc.variables['LON'][:]
latdrifter3 = nc.variables['LAT'][:]
nc.close()
nc=netcdf.Dataset(drifterfile4)
londrifter4 = nc.variables['LON'][:]
latdrifter4 = nc.variables['LAT'][:]
nc.close()
nc=netcdf.Dataset(drifterfile5)
londrifter5 = nc.variables['LON'][:]
latdrifter5 = nc.variables['LAT'][:]
nc.close()
nc=netcdf.Dataset(drifterfile6)
londrifter6 = nc.variables['LON'][:]
latdrifter6 = nc.variables['LAT'][:]
nc.close()

nc=netcdf.Dataset(sstdir+sstfile)
sst = nc.variables['SST'][:].squeeze()
nc.close()

nc=netcdf.Dataset(sstdir+sstcoord)
lonsst = nc.variables['lon'][:]
latsst = nc.variables['lat'][:]
nc.close()

# ----------------------------------------------------------------------------
# Create the figure
fig=plt.figure() 
ax = plt.subplot(111)


contoursst=plt.pcolor(lonsst, latsst,sst,cmap=cmapsst,\
                          norm=normsst,alpha=1.0,zorder=3)

# Drifter trajectories 
plt.plot(londrifter1,latdrifter1,'.-',color='k',ms=2,lw=2,label='MLI drifter',zorder=4)
plt.plot(londrifter2,latdrifter2,'.-',color='r',ms=2,lw=2,label='MLI drifter',zorder=4)
plt.plot(londrifter3,latdrifter3,'.-',color='y',ms=2,lw=2,label='MLI drifter',zorder=4)
plt.plot(londrifter4,latdrifter4,'.-',color='b',ms=2,lw=2,label='MLI drifter',zorder=4)
plt.plot(londrifter5,latdrifter5,'.-',color='g',ms=2,lw=2,label='MLI drifter',zorder=4)
plt.plot(londrifter6,latdrifter6,'.-',color='c',ms=2,lw=2,label='MLI drifter',zorder=4)


# Overlay visible image
plt.imshow(arr[:3,:,:].transpose((1, 2, 0)), extent=extent,zorder=2,alpha=0.9)

# Coastline
plt.plot(lonc,latc,color='k',lw=0.5,zorder=4)

### Add text
##xt,yt=1.5,39.15
##plt.text(xt,yt,'Ibiza',ha='center',fontname='Times New Roman',fontsize=16,bbox=dict(fc='w', ec=None, lw=0.5,alpha=0.9, pad=10.0),zorder=4)
##xt,yt=1.55,38.545
##plt.text(xt,yt,'Formentera',ha='center',fontname='Times New Roman',fontsize=16,bbox=dict(fc='w', ec=None, lw=0.5,alpha=0.9, pad=10.0),zorder=4)
##plt.text(1.02,39.4,'SARAL/AltiKa track',ha='left',fontname='Times New Roman',fontsize=14,bbox=dict(fc='w', ec=None, lw=0.5,alpha=0.7, pad=10.0),zorder=4)
##plt.text(0.27,38.2,'Drifter',ha='left',fontname='Times New Roman',fontsize=14,bbox=dict(fc='w', ec=None, lw=0.5,alpha=0.7, pad=10.0),zorder=4)
##plt.text(0.9,39.05,'Radar velocities',ha='center',fontname='Times New Roman',fontsize=14,bbox=dict(fc='w', ec=None, lw=0.5,alpha=0.7, pad=10.0),zorder=4)
##plt.text(1.05,38.2,'Deep glider',ha='left',fontname='Times New Roman',fontsize=14,bbox=dict(fc='w', ec=None, lw=0.5,alpha=0.7, pad=10.0),zorder=4)



if doplotsst==1:
    cbar=fig.colorbar(contoursst,orientation='horizontal',\
              fraction=0.08,shrink=0.9,aspect=15,pad=0.02)
    cbar.set_label(r'($^{\circ}$C)',fontname='Times New Roman',fontsize=16,rotation=0)
    cbar.set_ticks(newticks)



ax.set_xlim(lonmin,lonmax)
ax.set_ylim(latmin,latmax)
ax.set_xticklabels('')
ax.set_yticklabels('')
ax.xaxis.tick_top()

    
plt.savefig(figdir+figname+figtype2, dpi=300, facecolor='w', edgecolor='w',
         transparent=False, bbox_inches='tight', pad_inches=0.1)
####plt.savefig(figdir+figname+figtype2, dpi=300, facecolor='w', edgecolor='w',
####         transparent=False, bbox_inches='tight', pad_inches=0.1)
#plt.show()
plt.close()
