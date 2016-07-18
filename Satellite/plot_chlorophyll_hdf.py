#!/usr/bin/env python
#
# plot_chlorophyll_hdf.py
# 
# Plot the chlorophyll concentration
# 
# ctroupin, December 2013
# data downloaded from
# http://oceancolor.gsfc.nasa.gov
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


datadir='/home/ctroupin/DataOceano/Chlorophyll/'
figdir=datadir+'figures/'

databasefile = 'A*.L2_LAC_OC'
figname='chloro_'
figtype='.png'

basemap_resolution = 'l'
deg2km = 111125.1               

#------------------------------------------------------------------------------------
cmap=plt.cm.hot_r
cmap2=plt.cm.RdBu_r

vmin = 0.
vmax = 15.
vmin2=-1.25
vmax2=1.25
dvar=2.5
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

### canary
##lonmin = -13
##lonmax = -9
##latmin = 29
##latmax = 33

# West Med Sea
lonmin=-1.0
lonmax=5.0
latmin=38.0
latmax=41.0

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
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

#------------------------------------------------------------------------------------

# prepare time
deltatime=(dt.datestr2num('1990,1,1')-dt.datestr2num('1970,1,1'))*86400
filelist=sorted(glob.glob(datadir+databasefile))



# Loop on the files
for datafile in filelist:

    figbasename = os.path.basename(datafile)
    figbasename = figbasename[:-3]    
    
    print figbasename



##    #------------------------------------------------------------------------------------
##    # Select sub-region
##    goodlon = np.nonzero(np.logical_and(lon<=lonmax+dlon,lon>=lonmin-dlon))
##    goodlon = goodlon[0]
##    if goodlon.size !=0:
##        lat = lat[goodlon]
##        lon = lon[goodlon]
##        windspeed = windspeed[goodlon]
##        winddirection = -winddirection[goodlon]+90
##        goodlat = np.nonzero(np.logical_and(lat<=latmax+dlat,lat>=latmin-dlat))
##        goodlat = goodlat[0]
##        if goodlat.size !=0:
##            lat = lat[goodlat]
##            lon = lon[goodlat]
##            windspeed = windspeed[goodlat]
##            winddirection = winddirection[goodlat]
##
##            #------------------------------------------------------------------------------------
##            # Plot velocity field
##            if doplotvel ==1:
##
##                uwind=windspeed*np.cos(np.deg2rad(winddirection))
##                vwind=windspeed*np.sin(np.deg2rad(winddirection))
##                uwind=np.ma.masked_where(uwind==uwind.fill_value,uwind)
##                vwind=np.ma.masked_where(vwind==vwind.fill_value,vwind)
##                uwind.data[uwind.data==uwind.data.min()]=0
##                vwind.data[vwind.data==vwind.data.min()]=0
##                windnorm=np.sqrt(uwind*uwind+vwind*vwind)
##                x,y=m(lon,lat)
##                # Make the plot
##                fig=plt.figure()
##                ax = fig.add_subplot(111)
##                m.ax=ax
##
##                #contour=m.contourf(x,y,windspeed,levels2plot,cmap=cmap,norm=norm,extend='max')
##
##                Q=m.quiver(x,y,uwind,vwind,windnorm,\
##                           units='width',scale=250,width=0.003,norm=norm,cmap=cmap)
##                k = plt.quiverkey(Q, .9, 1.05, 10, r'$10\, ms^{-1}$',labelpos='W',
##                               fontproperties={'weight':'bold','size':'16'},color='k')
##                
##                # Add grid, coastline and continent
##                m.drawcoastlines(ax=ax)
##                m.fillcontinents(color='.15', ax=ax)
##                
##                meridians=np.arange(lonmin,lonmax,dlon)
##                parallels=np.arange(latmin,latmax,dlat)
##                m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,\
##                                fontname='Times New Roman')
##                m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,\
##                                fontname='Times New Roman')
##
##                # Add the colorbar
##                cbar=fig.colorbar(Q,cmap=cmap,orientation='vertical',
##                                  pad=0.025,aspect=12,shrink=0.65,norm=norm,extend='max')
##                cbar.set_label(r'$ms^{-1}$',fontname='Times New Roman',fontsize=18)
##                cbar.set_ticks(bounds)
##                cbar.ax.set_yticklabels(bounds,fontname='Times New Roman',fontsize=16)
##
##                plt.title(figtitle,fontname='Times New Roman',fontsize=24)
##                ## Export figure and display it
##                plt.savefig(figdir1+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
##                         transparent=False, bbox_inches='tight', pad_inches=0.1)
##                plt.show()
##                plt.close()
##        else:
##            print 'no data in your area'
##            #------------------------------------------------------------------------------------

