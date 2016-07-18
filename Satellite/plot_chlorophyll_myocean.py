#!/usr/bin/env python
#
# plot_chlorophyll_myocean.py
# 
# Plot the wind obtained from MyOcean
# 
# ctroupin, November 2013
# data downloaded from
# ftp://myocean.artov.isac.cnr.it/Core
# ---------------------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
from matplotlib.colors import LogNorm
import time
from mpl_toolkits.basemap import Basemap

doplot=1

chlorodir = '/home/ctroupin/DataOceano/Chlorophyll/MyOcean/'
figdir=chlorodir+'figures/'

chlorobasefile = '201402*.nc'
figtype='.png'

basemap_resolution = 'h'
deg2km = 111125.1               

#------------------------------------------------------------------------------------
cmap=plt.cm.jet

vmin = 0.
vmax = 15.
dvar=2.5

bounds = 10.**np.arange(-2.,1.,1.)
levels2plot = 10.**np.arange(-2.,1.,.05)
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
deltatime=(dt.datestr2num('1981,1,1')-dt.datestr2num('1970,1,1'))*86400


filelist=sorted(glob.glob(chlorodir+chlorobasefile))

# Loop on the files
for chlorofile in filelist:

    figbasename = os.path.basename(chlorofile)
    figbasename = figbasename[:-3]    
    
    # Load chlorophyll concentration
    nc=netcdf.Dataset(chlorofile)
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    chloro = nc.variables['CHL'][:]
    chlorotime = nc.variables['time'][:]
    nc.close()

    chloro=np.squeeze(chloro)
    lon[lon>180]=lon[lon>180]-360.0

    lon=lon.flatten()
    lat=lat.flatten()

    figtime=time.gmtime(chlorotime[0]+deltatime)
    figtitle=str(figtime.tm_year)+'-'+str(figtime.tm_mon)+'-'+str(figtime.tm_mday)

    # Select sub-region
    goodlon = np.nonzero(np.logical_and(lon<=lonmax+dlon,lon>=lonmin-dlon))
    goodlat = np.nonzero(np.logical_and(lat<=latmax+dlat,lat>=latmin-dlat))
    goodlat = goodlat[0][:]
    goodlon = goodlon[0][:]

    lat = lat[goodlat]
    lon = lon[goodlon]
    
    chloro=chloro[goodlat,:]
    chloro=chloro[:,goodlon]

    # Mask
    np.ma.masked_where(chloro==chloro.fill_value,chloro)

    llon,llat=np.meshgrid(lon,lat)
    
    #------------------------------------------------------------------------------------
    # Plot chlorophyll field
    if doplot ==1:
        
        # Make the plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        m.ax=ax
        x,y=m(llon,llat)
        contour=m.contourf(x,y,chloro,levels2plot,cmap=cmap,norm = LogNorm())

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
        cbar=fig.colorbar(contour,cmap=cmap,orientation='vertical',
                          pad=0.025,aspect=12,shrink=0.65,extend='max')
        cbar.set_label(r'$mg\,m^{-3}$',fontname='Times New Roman',fontsize=18)
        cbar.set_ticks(bounds)
        cbar.ax.set_yticklabels(bounds,fontname='Times New Roman',fontsize=16)

        plt.title(figtitle,fontname='Times New Roman',fontsize=24)
        ## Export figure and display it
        plt.savefig(figdir+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        plt.close()





