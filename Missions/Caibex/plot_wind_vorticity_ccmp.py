#!/usr/bin/env python
#
# 
# Compute and plot the wind vorticity
# obtained from CCMP
#
# ctroupin, November 2013
# --------------------------------------------------------------


import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
import time
from mpl_toolkits.basemap import Basemap


doplotvel=0
doplotcurl=1

winddir = '/home/ctroupin/DataOceano/Wind/CCMP/data/'
figdir1='/home/ctroupin/DataOceano/Wind/CCMP/figures/wind/'
figdir2='/home/ctroupin/DataOceano/Wind/CCMP/figures/windcurl3/'

windbasefile = 'analysis_*_v11l30flk.nc'

figname='windcurl_'
figtype='.png'

filelist = sorted(glob.glob(winddir+windbasefile))

basemap_resolution = 'h'
deg2km = 111125.1               

cmap=plt.cm.jet
cmap2=plt.cm.RdBu_r

vmin = 0.
vmax = 10.
vmin2=-1.25
vmax2=1.25
dvar=2.
dvar2=0.25
bounds = np.arange(vmin,vmax+0.01,dvar)
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.01)
levels2plot = np.arange(vmin,vmax+0.001,dvar/10)

bounds2 = np.arange(vmin2,vmax2+0.01,dvar2)
bounds2[len(bounds2)/2] = 0
norm2 = colors.Normalize(vmin=vmin2,vmax=vmax2+0.01)
levels2plot2 = np.arange(vmin2,vmax2+0.001,dvar2/10)

#------------------------------------------------------------------------------------
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
# Create figure directories if necessary
if not(os.path.exists(figdir1)):
    os.makedirs(figdir1)
if not(os.path.exists(figdir2)):
    os.makedirs(figdir2)
#------------------------------------------------------------------------------------
# Load wind

# Load coordinates from 1st file
nc=netcdf.Dataset(filelist[0])
lon = nc.variables['lon'][:]-360
lat = nc.variables['lat'][:]
nc.close()

#------------------------------------------------------------------------------------
# Select sub-region
goodlon = np.nonzero(np.logical_and(lon<=lonmax+dlon,lon>=lonmin-dlon))
goodlat = np.nonzero(np.logical_and(lat<=latmax+dlat,lat>=latmin-dlat))
goodlat = goodlat[0][:]
goodlon = goodlon[0][:]

lat = lat[goodlat]
lon = lon[goodlon]

llon,llat=np.meshgrid(lon,lat)
x,y = m(llon,llat)

#------------------------------------------------------------------------------------        
# Compute metrics and Coriolis frq
lon,lat=np.meshgrid(lon,lat)
dlonx,dlony=np.gradient(lon)
dlatx,dlaty=np.gradient(lat)
dx=dlony*np.cos(np.deg2rad(lat))*deg2km
dy=dlatx*deg2km
f=4*np.pi*np.sin(np.deg2rad(lat))/86400

# Loop on the files
for windfile in filelist:

    figbasename = os.path.basename(windfile)
    figbasename = figbasename[:-3]    
    figbasename2 = figbasename+'_curl'
    
    # Load wind
    nc=netcdf.Dataset(windfile)
    uwind = nc.variables['uwnd'][:]
    vwind = nc.variables['vwnd'][:]
    windtime = nc.variables['time'][:]
    nc.close()

    # Average to get daily fields
    uwind=np.mean(uwind,axis=0)
    vwind=np.mean(vwind,axis=0)

    uwind = uwind[goodlat]
    uwind = uwind[:,goodlon]
    vwind = vwind[goodlat]
    vwind = vwind[:,goodlon]

    # Construct year-month-day
    deltatime=dt.datestr2num('1987,1,1')-dt.datestr2num('1970,1,1')
    ttime=time.gmtime((deltatime+(windtime[0]/24))*86400)
    figtitle=str(ttime.tm_year)+'-'+str(ttime.tm_mon)+'-'+str(ttime.tm_mday)

    #------------------------------------------------------------------------------------
    # Plot velocity field
    if doplotvel ==1:
    
        uvnorm= np.sqrt(uwind*uwind+vwind*vwind)

        # Make the plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        m.ax=ax

        contour=m.contourf(x,y,uvnorm,levels2plot,cmap=cmap,norm=norm,extend='max')
        Q=m.quiver(x,y,uwind,vwind,\
                   width=0.003,color='k')
        k = plt.quiverkey(Q, .95, 1.075, 10, r'$10\, ms^{-1}$',labelpos='S',
                       fontproperties={'weight':'bold','size':'16'},color='k')
        
        # Add grid, coastline and continent
        m.drawcoastlines(ax=ax)
        m.fillcontinents(color='.15', ax=ax)
        
        meridians=np.arange(lonmin,lonmax,dlon)
        parallels=np.arange(latmin,latmax,dlat)
        m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,\
                        fontname='Times New Roman')
        m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,\
                        fontname='Times New Roman')

        # Add the colorbar
        cbar=fig.colorbar(contour,cmap=cmap,orientation='vertical',pad=0.025,aspect=15)
        cbar.set_label(r'$ms^{-1}$',fontname='Times New Roman',fontsize=18)
        cbar.set_ticks(bounds)
        cbar.ax.set_yticklabels(bounds,fontname='Times New Roman',fontsize=16)

        plt.title(figtitle,fontname='Times New Roman',fontsize=24)
        ## Export figure and display it
        plt.savefig(figdir1+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        plt.close()

    #------------------------------------------------------------------------------------
    # Plot wind curl
    if doplotcurl ==1:
    
        
        # Compute wind curl
        dux,duy=np.gradient(uwind)
        dvx,dvy=np.gradient(vwind)
        vort=(dvy/dx-dux/dy)/f

        
        # Make the plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        m.ax=ax

        contour2=m.contourf(x,y,vort,levels2plot2,cmap=cmap2,norm=norm2,extend='both')
        Q=m.quiver(x,y,uwind,vwind,\
                   width=0.003,color='k')
        k = plt.quiverkey(Q, .95, 1.075, 10, r'$10\, ms^{-1}$',labelpos='S',
                       fontproperties={'weight':'bold','size':'16'},color='k')
        
        # Add grid, coastline and continent
        m.drawcoastlines(ax=ax)
        m.fillcontinents(color='.15', ax=ax)
        
        meridians=np.arange(lonmin,lonmax,dlon)
        parallels=np.arange(latmin,latmax,dlat)
        m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,\
                        fontname='Times New Roman')
        m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,\
                        fontname='Times New Roman')

        # Add the colorbar

        cbar=fig.colorbar(contour2,cmap=cmap2,orientation='vertical',pad=0.025,aspect=15)
        cbar.set_label(r'$\xi/f$',fontname='Times New Roman',fontsize=18)
        cbar.set_ticks(bounds2)
        cbar.ax.set_yticklabels(bounds2,fontname='Times New Roman',fontsize=16)

        plt.title(figtitle,fontname='Times New Roman',fontsize=24)
        ## Export figure and display it
        plt.savefig(figdir2+figbasename2+figtype, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        plt.close()
    
