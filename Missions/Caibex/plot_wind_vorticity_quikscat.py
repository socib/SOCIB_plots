#!/usr/bin/env python
#
# 
# Compute and plot the wind vorticity
# obtained from QuikSCAT
#
# ctroupin, November 2013

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
import time, datetime
from mpl_toolkits.basemap import Basemap
import scipy.io

winddir = '/data_local/Satellite/Wind/QuikScat/Results/'
windfile = 'qscat4dineof_300e360e0n60n_filtered.nc'


figdir='/data_local/Satellite/Wind/QuikScat/figures/'
figname='windcurl_'
figtype='.png'


cmap=plt.cm.RdBu_r
basemap_resolution = 'l'
deg2km = 111125.1               

# Region of interest
lonmin = -13
lonmax = -9
latmin = 29
latmax = 33
dlon=1
dlat=1

m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                lat_ts=0.5*(latmin+latmax),\
                resolution=basemap_resolution)

meridians=np.arange(lonmin,lonmax,dlon)
parallels=np.arange(latmin,latmax,dlat)

#------------------------------------------------------------------------------------

# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)
#------------------------------------------------------------------------------------
# Load wind

##mat = scipy.io.loadmat(winddir+windfile)
##lon = mat['llon_wind'][:]
##lat = mat['llat_wind'][:]
##u = mat['uwindmean'][:]
##v = mat['vwindmean'][:]

# Load coordinates from 1st file
nc=netcdf.Dataset(winddir+windfile)
lon = nc.variables['lon'][:]-360
lat = nc.variables['lat'][:]
uwind = nc.variables['uwind'][:]
vwind = nc.variables['vwind'][:]
time = nc.variables['time'][:]
nc.close()

# Select period of interest
firstday = dt.datestr2num('2010,8,15')
lastday = dt.datestr2num('2010,9,6')
goodtime = np.nonzero(np.logical_and(time>=firstday,time<=lastday))
goodtime=goodtime[0]

uwind=uwind[goodtime,:,:]
vwind=vwind[goodtime,:,:]
time2=time[goodtime]

# Mask 
uwind = np.ma.masked_array(uwind,uwind==uwind.fill_value)
vwind = np.ma.masked_array(vwind,uwind==vwind.fill_value)

#------------------------------------------------------------------------------------                    
# Compute vorticity
lon,lat=np.meshgrid(lon,lat)
dlonx,dlony=np.gradient(lon)
dlatx,dlaty=np.gradient(lat)
dx=dlony*np.cos(np.deg2rad(lat))*deg2km
dy=dlatx*deg2km

# Compute Coriolis frequency
f=4*np.pi*np.sin(np.deg2rad(lat))/86400

# Loop on time
for tt in range(0,1):#len(time)-1):
    dux,duy=np.gradient(np.squeeze(uwind[tt,:,:]))
    dvx,dvy=np.gradient(np.squeeze(vwind[tt,:,:]))
    vort=(dvy/dx-dux/dy)/f

    #------------------------------------------------------------------------------------
    # Make the plot

    fig=plt.figure()
    ax = fig.add_subplot(111)
    lon2,lat2=m(lon,lat)
    Q=m.quiver(lon2,lat2,np.squeeze(uwind[tt,:,:]),np.squeeze(vwind[tt,:,:]),\
                     units='width',scale=200,width=0.003,zorder=3)

    k = plt.quiverkey(Q,1.0,1.06, 5, r'$5\, ms^{-1}$',labelpos='S',
               fontproperties={'weight':'bold','size':'18'},color='k')
    
    pcm=m.pcolormesh(lon2,lat2,vort,zorder=2,cmap=cmap)
    cbar=fig.colorbar(pcm)
    cbar.drawedges

    m.drawmeridians(meridians,labels=[0, 0, 1, 0],linewidth=0.5,fontsize=18,\
                    fontname='Times New Roman')
    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,\
                    fontname='Times New Roman')
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='.15', ax=ax)

##    plt.savefig(figdir+figname+str(tt)+figtype, dpi=300, facecolor='w', edgecolor='w',
##                 transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close()
