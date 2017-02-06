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
datadir1='/home/ctroupin/DataOceano/CAIBEX_campaign/data/CTD/data4diva/'
datadir2='/home/ctroupin/DataOceano/CAIBEX_campaign/data/Seasor/data4diva/'
resultdir='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/results/Horizontal/V2/'
resultbasefile=varname+'1001*.nc'
figdir='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/figures/results/horizontal/V2/'
resultfilelist=sorted(glob.glob(resultdir+resultbasefile))
contourfile='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/contours/contour.depth'


#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'f'
maxerror=0.25    # Maximal relative error for the plot


# region of interest
lonmin = -12.25 
lonmax = -9.5 
latmin = 29.75
latmax = 31.75

levels2plot=np.arange(12,23,0.25)

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

depthlist=np.loadtxt(contourfile)

i=0
# Loop on result files
for resultfile in resultfilelist:

    # Figure name
    figbasename = os.path.basename(resultfile).split('.')[0]
    print figbasename
    num=figbasename[-5:]

    datafile1=datadir1+varname+'_CTD_'+num+'.dat'
    datafile2=datadir2+varname+'_seasoar_'+num+'.dat'
    
    londata1,latdata1=np.loadtxt(datafile1,usecols=(0,1),unpack=True)
    londata2,latdata2=np.loadtxt(datafile2,usecols=(0,1),unpack=True)
    
    # Load results from file
    nc=netcdf.Dataset(resultfile)
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
    x2,y2 = m(londata2,latdata2)

    contour=m.contour(x,y,field,levels2plot,colors='k',zorder=2)
    plt.clabel(contour,levels2plot[0::2],inline=1,fmt='%1.1f',fontsize=14,fontname='Times New Roman')

    plt.plot(x1,y1,'ko',ms=3)
    plt.plot(x2,y2,'ks',ms=3)
    
    # Add grid, coastline and continent
    m.drawcoastlines(ax=ax,zorder=4)
    m.fillcontinents(color='0.7', ax=ax,zorder=3)
    m.drawmeridians(lon2plot,labels=[0, 0, 0, 1],fontsize=18,fontname='Times New Roman')
    m.drawparallels(lat2plot,labels=[1, 0, 0, 0],fontsize=18,fontname='Times New Roman')

    titletext = str(int(depthlist[int(num[-2:])])) + ' m'
    i+=1
    plt.title(titletext,fontname='Times New Roman',fontsize=20)

    # Export figure and display it
    plt.savefig(figdir+figbasename, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)

   # plt.show()
    plt.close()

