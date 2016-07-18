#!/usr/bin/env python
#
# 
# Plot the L2 wind fields obtained from MetOp-a or -b
# 
# ctroupin, November 2013
# data downloaded from
# ftp://podaac-ftp.jpl.nasa.gov/allData/ascat/preview/L2/metop_a/coastal_opt/
# ---------------------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
import time
from mpl_toolkits.basemap import Basemap


windbasedir = '/home/ctroupin/DataOceano/Wind/ASCAT/data/2013/'
figdir1='/home/ctroupin/DataOceano/Wind/ASCAT/figures/'
windbasefile = 'ascat_*.nc'

# Days of the year 
firstday,lastday=213,214

figname1='windvec_'
figtype='.png'

basemap_resolution = 'h'
deg2km = 111125.1               

#------------------------------------------------------------------------------------

# colormap
cmap=plt.cm.hot_r

# Wind speed limits
vmin = 0.
vmax = 15.
# Spacing between labels on the colormap
dvar = 2.5

bounds = np.arange(vmin,vmax+0.01,dvar)
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.01)

#------------------------------------------------------------------------------------
# Region of interest
# West Med Sea
lonmin=-1.0
lonmax=5.0
latmin=38.0
latmax=41.0

dlon=1
dlat=1

# Prepare projection
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                lat_ts=0.5*(latmin+latmax),\
                resolution=basemap_resolution)

# Meridians and parallels to plot
meridians=np.arange(lonmin,lonmax,dlon)
parallels=np.arange(latmin,latmax,dlat)

#------------------------------------------------------------------------------------
# Create figure directories if necessary
if not(os.path.exists(figdir1)):
    os.makedirs(figdir1)
#------------------------------------------------------------------------------------

# prepare time (reference to January 1st, 1990)
deltatime=(dt.datestr2num('1990,1,1')-dt.datestr2num('1970,1,1'))*86400


# Loop on the selected days
for dday in range(firstday,lastday,):
    winddir=windbasedir+str(dday)+'/'
    filelist=sorted(glob.glob(winddir+windbasefile))

    # Loop on the files
    for windfile in filelist:

        figbasename = os.path.basename(windfile)
        figbasename = figbasename[:-3]    
        figbasename2 = figbasename+'_curl'
        
        # Load wind
        nc=netcdf.Dataset(windfile)
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        windspeed = nc.variables['wind_speed'][:]
        winddirection = nc.variables['wind_dir'][:]
        windtime = nc.variables['time'][:]
        nc.close()

        # Change latitudes to have them between -180 and 180
        lon[lon>180]=lon[lon>180]-360.0

        # Put into vector format
        lon=lon.flatten()
        lat=lat.flatten()
        windspeed=windspeed.flatten()
        winddirection=winddirection.flatten()

        # Compute time 
        figtime=time.gmtime(windtime[0,0]+deltatime)
        figtitle=str(figtime.tm_year)+'-'+str(figtime.tm_mon)+'-'+str(figtime.tm_mday)
        
        #------------------------------------------------------------------------------------
        # Select sub-region and check if data fall inside

        # Longitudes
        goodlon = np.nonzero(np.logical_and(lon<=lonmax+dlon,lon>=lonmin-dlon))
        goodlon = goodlon[0]
        if goodlon.size !=0:
            lat = lat[goodlon]
            lon = lon[goodlon]
            windspeed = windspeed[goodlon]

            # Correction because wind in files is defined with 0 degree northward.
            winddirection = -winddirection[goodlon]+90

            # Latitude
            goodlat = np.nonzero(np.logical_and(lat<=latmax+dlat,lat>=latmin-dlat))
            goodlat = goodlat[0]

            # Check if data exist
            if goodlat.size !=0:
                lat = lat[goodlat]
                lon = lon[goodlat]
                windspeed = windspeed[goodlat]
                winddirection = winddirection[goodlat]

                #------------------------------------------------------------------------------------
                # Compute components using wind direction
                uwind=windspeed*np.cos(np.deg2rad(winddirection))
                vwind=windspeed*np.sin(np.deg2rad(winddirection))

                # Mask baed values
                uwind=np.ma.masked_where(uwind==uwind.fill_value,uwind)
                vwind=np.ma.masked_where(vwind==vwind.fill_value,vwind)

                # Put minimum to zero (should not be necessary; to check)
                uwind.data[uwind.data==uwind.data.min()]=0
                vwind.data[vwind.data==vwind.data.min()]=0
                windnorm=np.sqrt(uwind*uwind+vwind*vwind)
                x,y=m(lon,lat)
                
                # Make the plot
                fig=plt.figure()
                ax = fig.add_subplot(111)
                m.ax=ax

                # Vectors
                Q=m.quiver(x,y,uwind,vwind,windnorm,units='width',scale=250,width=0.003,norm=norm,cmap=cmap)

                # Reference vector
                k = plt.quiverkey(Q, .9, 1.05, 10, r'$10\, ms^{-1}$',labelpos='W',
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
                cbar=fig.colorbar(Q,cmap=cmap,orientation='vertical',
                                  pad=0.025,aspect=12,shrink=0.65,norm=norm,extend='max')
                cbar.set_label(r'$ms^{-1}$',fontname='Times New Roman',fontsize=18)
                cbar.set_ticks(bounds)
                cbar.ax.set_yticklabels(bounds,fontname='Times New Roman',fontsize=16)

                plt.title(figtitle,fontname='Times New Roman',fontsize=24)
                
                ## Export figure and display it
                plt.savefig(figdir1+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
                         transparent=False, bbox_inches='tight', pad_inches=0.1)

                # Uncomment next line once everything is ok
                #plt.show()
                plt.close()
            else:
                print 'no data in your area'
                #------------------------------------------------------------------------------------
                # 



