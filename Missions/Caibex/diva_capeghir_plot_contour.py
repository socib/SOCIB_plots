#!/usr/bin/env python
# diva_plot_contour.py
#------------------------------------------------------------------------------

import os, glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import linecache

# Clean 
os.system('clear')

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'l'

# File and directory names
contourdir='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/contours/'
figdir='/home/ctroupin/DataOceano/CAIBEX_campaign/Diva/figures/contours/'
contourfiles= sorted(glob.glob(contourdir+'coast.cont*'))


# Figure name
figbasename = 'contour_capeghir'

# Region of interest        
lonmin=-13.
lonmax=-9.
latmin=29
latmax=32

# Spacing between x/y labels
dlon = 1
dlat = 1

# Unit of the variable to plot
unitname = ' '

#------------------------------------------------------------------------------------

# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

i=0
for file2load in contourfiles:
    print file2load
    i+=1
    figname=figdir+figbasename+str(i).zfill(3)
    print figname

                   
    
    # Count number of files and contours
    ncontours=int(linecache.getline(file2load,1))
    nlines = sum(1 for line in open(file2load))

    # Prepare figure and projection
    fig=plt.figure()
    ax = fig.add_subplot(111)

    m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
                    urcrnrlon=lonmax,urcrnrlat=latmax,  \
                    lat_ts=0.5*(lonmin+lonmax),\
                    resolution=basemap_resolution)
    m.ax=ax

    # Initialise line to read number
    linenum=2

    # Start loop on the contours
    for n in (np.arange(0,ncontours)):

        # Number of points in the current contour
        npoints=int(linecache.getline(file2load,linenum))
        
        nskiplines = linenum+npoints
        # Load coordinates (npoints lines to be read)
        coord=np.genfromtxt(file2load,skip_header=linenum,skip_footer=nlines-nskiplines)
        lon=coord[:,0]
        lat=coord[:,1]

        x,y = m(lon, lat)
        m.plot(x,y,'k',lw=1)
        
        # Update line number
        # (taking into accountnumber of points read)
        linenum=nskiplines+1

    # Add grid, coastline and continent
    #m.fillcontinents(color='black', ax=ax)
    meridians=np.arange(lonmin,lonmax,dlon)
    parallels=np.arange(latmin,latmax,dlat)
    m.drawmeridians(meridians,labels=[0, 0, 1, 0],fontsize=18,fontname='Times New Roman')
    m.drawparallels(parallels,labels=[1, 0, 0, 0],fontsize=18,fontname='Times New Roman')

    # Export figure and display it
    plt.savefig(figname, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)

    #plt.show()
    plt.close()
