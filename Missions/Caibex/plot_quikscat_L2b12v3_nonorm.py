#!/usr/bin/env python
#
# 
# Plot the wind obtained from QuikSCAT
# 
# ctroupin, November 2013
# data downloaded from ftp://podaac.jpl.nasa.gov/allData/quikscat/L2B12/v3/

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
import time
import gmtColormap
from mpl_toolkits.basemap import Basemap
from scipy import interpolate

doplotvel=0
doplotcurl=1

windbasedir = '/home/ctroupin/DataOceano/Satellite/QuikSCAT/netcdf/L2B12v3/data/'
figdir1='/home/ctroupin/DataOceano/Satellite/QuikSCAT/netcdf/L2B12v3/figures/windcolor3/'
figdir2='/home/ctroupin/DataOceano/Satellite/QuikSCAT/netcdf/L2B12v3/figures/windcurl_no_normalized/'

windbasefile = 'qs_l2b*.nc'

firstday,lastday=237,246
#firstday,lastday=240,241

figname='windcurl_'
figtype='.png'

basemap_resolution = 'f'
deg2km = 111125.1               

#------------------------------------------------------------------------------------
##cmap=plt.cm.autumn_r
##cmap2=plt.cm.RdBu_r

cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps')
cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)

cdict= gmtColormap.gmtColormap('temp_19lev','/home/ctroupin/Software/Python/GMT_colormaps')
cmap2 = colors.LinearSegmentedColormap('my_colormap',cdict,256)

vmin = 2.
vmax = 10.
vmin2=-1.5
vmax2=1.5
dvar=2.
dvar2=0.5

bounds = np.arange(vmin,vmax+0.01,dvar)
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.01)
levels2plot = np.arange(vmin,vmax+0.001,dvar/10)

bounds2 = np.arange(vmin2,vmax2+0.01*dvar2,dvar2)
bounds2[len(bounds2)/2] = 0
norm2 = colors.Normalize(vmin=vmin2,vmax=vmax2+0.01*dvar2)
levels2plot2 = np.arange(vmin2,vmax2+0.01*dvar2,dvar2/10)

#------------------------------------------------------------------------------------
# Region of interest
lonmin = -13
lonmax = -8.9
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

# Create regular grid on which wind will be interpolated
# ------------------------------------------------------

dloninterp=0.10
dlatinterp=0.10
lon2interp=np.arange(lonmin-2*dloninterp,lonmax+2*dloninterp,dloninterp)
lat2interp=np.arange(latmin-2*dlatinterp,latmax+2*dlatinterp,dlatinterp)
llon2interp,llat2interp=np.meshgrid(lon2interp,lat2interp)

# Compute metrics and Coriolis frq
dlonx,dlony=np.gradient(llon2interp)
dlatx,dlaty=np.gradient(llat2interp)
dx=dlony*np.cos(np.deg2rad(llat2interp))*deg2km
dy=dlatx*deg2km
f=4*np.pi*np.sin(np.deg2rad(llat2interp))/86400

#------------------------------------------------------------------------------------
# Create figure directories if necessary
if not(os.path.exists(figdir1)):
    os.makedirs(figdir1)
if not(os.path.exists(figdir2)):
    os.makedirs(figdir2)
#------------------------------------------------------------------------------------

# prepare time
deltatime=(dt.datestr2num('1999,1,1')-dt.datestr2num('1970,1,1'))*86400


# Loop on the days
for dday in range(firstday,lastday,):
    winddir=windbasedir+str(dday)+'/'
    filelist=sorted(glob.glob(winddir+windbasefile))


    
    ###------------------------------------------------------------------------------------        
    ### Compute metrics and Coriolis frq
    ##lon,lat=np.meshgrid(lon,lat)
    ##dlonx,dlony=np.gradient(lon)
    ##dlatx,dlaty=np.gradient(lat)
    ##dx=dlony*np.cos(np.deg2rad(lat))*deg2km
    ##dy=dlatx*deg2km
    ##f=4*np.pi*np.sin(np.deg2rad(lat))/86400

    # Loop on the files
    for windfile in filelist:

        figbasename = os.path.basename(windfile)
        figbasename = figbasename[:-3]    
        figbasename2 = figbasename+'_curl'
        
        # Load wind
        nc=netcdf.Dataset(windfile)
        lon = nc.variables['lon'][:]-360
        lat = nc.variables['lat'][:]
        windspeed = nc.variables['retrieved_wind_speed'][:]
        winddirection = nc.variables['retrieved_wind_direction'][:]
        windtime = nc.variables['time'][:]
        nc.close()

        lon=lon.flatten()
        lat=lat.flatten()
        windspeed=windspeed.flatten()
        winddirection=winddirection.flatten()

        figtime=time.gmtime(windtime[0]+deltatime)
        figtitle=str(figtime.tm_year)+'-'+str(figtime.tm_mon)+'-'+str(figtime.tm_mday)
        #------------------------------------------------------------------------------------
        # Select sub-region
        goodlon = np.nonzero(np.logical_and(lon<=lonmax+dlon,lon>=lonmin-dlon))
        goodlon = goodlon[0]
        if goodlon.size !=0:
            lat = lat[goodlon]
            lon = lon[goodlon]
            windspeed = windspeed[goodlon]
            winddirection = -winddirection[goodlon]+90
            goodlat = np.nonzero(np.logical_and(lat<=latmax+dlat,lat>=latmin-dlat))
            goodlat = goodlat[0]
            if goodlat.size !=0:
                lat = lat[goodlat]
                lon = lon[goodlat]
                windspeed = windspeed[goodlat]
                winddirection = winddirection[goodlat]

                uwind=windspeed*np.cos(np.deg2rad(winddirection))
                vwind=windspeed*np.sin(np.deg2rad(winddirection))
                uwind=np.ma.masked_where(uwind==uwind.fill_value,uwind)
                vwind=np.ma.masked_where(vwind==vwind.fill_value,vwind)
                
                uwind.data[uwind.data==uwind.data.min()]=np.nan
                vwind.data[vwind.data==vwind.data.min()]=np.nan
                uwind=np.ma.masked_where(uwind.data==0.0,uwind)
                vwind=np.ma.masked_where(vwind.data==0.0,vwind)
                windnorm=np.sqrt(uwind*uwind+vwind*vwind)

                #------------------------------------------------------------------------------------
                # Plot velocity field
                if doplotvel ==1:

                    
                    x,y=m(lon,lat)
                    # Make the plot
                    fig=plt.figure()
                    ax = fig.add_subplot(111)
                    m.ax=ax

                    #contour=m.contourf(x,y,windspeed,levels2plot,cmap=cmap,norm=norm,extend='max')

                    Q=m.quiver(x,y,uwind,vwind,windnorm,\
                               units='width',scale=150,width=0.003,norm=norm,cmap=cmap)
                    k = plt.quiverkey(Q, .9, 1.075, 10, r'$10\, ms^{-1}$',labelpos='S',
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
                    cbar=fig.colorbar(Q,orientation='vertical',pad=0.025,aspect=15,norm=norm,extend='both')
                    #cbar.set_label(r'$ms^{-1}$',fontname='Times New Roman',fontsize=18)
                    cbar.set_ticks(bounds)
                    cbar.ax.set_yticklabels(bounds,fontname='Times New Roman',fontsize=16)

                    plt.title(figtitle,fontname='Times New Roman',fontsize=24)
                    ## Export figure and display it
                    plt.savefig(figdir1+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
                             transparent=False, bbox_inches='tight', pad_inches=0.1)
                    #plt.show()
                    plt.close()

                if doplotcurl ==1:

                    x,y=m(llon2interp,llat2interp)

                    # Need to reinterpolate wind on a regular grid
                    uwind_interp = interpolate.griddata((lon,lat), uwind, (llon2interp,llat2interp), method='linear')
                    vwind_interp = interpolate.griddata((lon,lat), vwind, (llon2interp,llat2interp), method='linear')                    
                    
                    # Compute wind curl
                    dux,duy=np.gradient(uwind_interp)
                    dvx,dvy=np.gradient(vwind_interp)
                    vort=(dvy/dx-dux/dy)

                    vort=np.ma.masked_where(np.logical_or(np.isnan(uwind_interp),uwind_interp==0.0),vort)
             
                     
                    # Make the plot
                    fig=plt.figure()
                    ax = fig.add_subplot(111)
                    m.ax=ax
             
                    contour2=m.contourf(x,y,vort*1e4,levels2plot2,norm=norm2,cmap=cmap2,extend='both')
                    
##                    Q=m.quiver(x,y,uwind_interp,vwind_interp,units='width',scale=150,width=0.003,color='k')
##                               
##                    k = plt.quiverkey(Q, .95, 1.075, 10, r'$10\, ms^{-1}$',labelpos='S',
##                                    fontproperties={'weight':'bold','size':'16'},color='k')
                     
                    # Add grid, coastline and continent
                    m.drawcoastlines(ax=ax)
                    m.fillcontinents(color='.15', ax=ax)
                     
                    meridians=np.arange(lonmin,lonmax,dlon)
                    parallels=np.arange(latmin,latmax,dlat)
                    m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18)
                    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18)
             
                    # Add the colorbar
             
                    cbar=fig.colorbar(contour2,cmap=cmap2,orientation='vertical',pad=0.025,aspect=15)
                    xt,yt=m(-8.85,33.1)
                    plt.text(xt,yt,r'$( \times 10^{4} s^{-1})$',fontsize=18)

                    #cbar.set_label(r'$\xi/f$',fontname='Times New Roman',fontsize=18)

                    cbar.set_ticks(bounds2)
                    cbar.ax.set_yticklabels(bounds2,fontname='Times New Roman',fontsize=16)
             
                    plt.title(figtitle,fontname='Times New Roman',fontsize=24)
                    #   Export figure and display it
                    plt.savefig(figdir2+figbasename2+'_V2'+figtype, dpi=300, facecolor='w', edgecolor='w',
                         transparent=False, bbox_inches='tight', pad_inches=0.1)
                    #plt.show()
                    plt.close()
            else:
                print 'no data in your area'
                #------------------------------------------------------------------------------------
                # 






    ##    
         
