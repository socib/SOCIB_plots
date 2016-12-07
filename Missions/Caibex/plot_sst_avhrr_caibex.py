#!/usr/bin/env python

# plot_sst_avhrr_caibex.py
#
#------------------------------------------------------------------------------

import os
import glob 
import numpy as np
import math
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import gmtColormap
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors

# Clean 
os.system('clear')

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'h'

# File and directory names
datadir='/media/ctroupin/Iomega_HDD/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
databasename='2009*.nc'
figdir='/media/ctroupin/Iomega_HDD/DataOceano/Satellite/AVHRR/NAR/SSTonly//figures_13w9w29n33n_python/'
coordfile=datadir+'coord_NAR_SST.nc'

# Figure extension
figtype1='.png'
figtype2='.eps'

# Colormap
cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps')
cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)


# Compute min and max values 
vmin=17.0
vmax=24.0
norm = colors.Normalize(vmin=vmin,vmax=vmax)
nlevels=150
levels2plot=np.arange(vmin,vmax+0.0001,(vmax-vmin)/(nlevels-1))
newticks=np.arange(np.floor(vmin),np.ceil(vmax)+0.0001,1.0)

    
# Region of interest        
lonmin=-13.0
lonmax=-9.0
latmin=29.0
latmax=33.0

# Limits in the matrices (determined emperically)
imin,imax,jmin,jmax=450,800,1500,1800

# Spacing between x/y labels
dlon = 1
dlat = 1

# Unit of the variable to plot
unitname = ' '

plt.rc('font', serif='Times New Roman')
font = {'family':'serif','size':18}
plt.rc('font',**font)


m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                lat_ts=0.5*(latmin+latmax),\
                resolution=basemap_resolution)

#------------------------------------------------------------------------------------

# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)
    

# Loop on files
filelist = sorted(glob.glob(datadir+databasename))
nfiles = len(filelist)

# Load coordinates from coordinate file
nc=netcdf.Dataset(coordfile)
lon = nc.variables['lon'][imin:imax,jmin:jmax]
lat = nc.variables['lat'][imin:imax,jmin:jmax]
nc.close()

print lon.min()
print lon.max()
print lat.min()
print lat.max()




##
### Select MDT in the sub-region
##goodlon = np.nonzero(np.logical_and(lon<=lonmax+dlon,lon>=lonmin-dlon))
##goodlat = np.nonzero(np.logical_and(lat<=latmax+dlat,lat>=latmin-dlat))
##goodlat = goodlat[0][:]
##goodlon = goodlon[0][:]
##
##lat = lat[goodlat]
##lat = lat[:,goodlon]
##lon = lon[goodlat]
##lon = lon[:,goodlon]
##
for SSTfile in filelist:
    filename= os.path.basename(SSTfile)
    figname=filename.split(".")[0]
    figtitle=str(figname[0:4])+'-'+str(figname[4:6])+'-'+str(figname[6:8])
                                                    
    print figname

    # Load data from current file
    nc=netcdf.Dataset(SSTfile)
    field = nc.variables['sst'][0,imin:imax,jmin:jmax]
    nc.close()

    # Select sub-region
    field = field.squeeze()
   

    # Apply transformation
   # field=1.0*(field.squeeze()-273.15)

    # Mask land values
    field = np.ma.masked_array(field,field==field.fill_value)

    newticks=newticks.round(2)

    # Make the plot
    fig=plt.figure()
    ax = fig.add_subplot(111)


    m.ax=ax

    
    x,y = m(lon,lat)

    pcm=m.pcolormesh(x,y,field,cmap=cmap,norm=norm)

    # Add grid, coastline and continent
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='.15', ax=ax)
    #m.bluemarble()

    meridians=np.arange(lonmin,lonmax,dlon)
    parallels=np.arange(latmin,latmax,dlat)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,\
                    fontname='Times New Roman')
    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,\
                    fontname='Times New Roman')

    # Add the colorbar
    cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
    cbar.set_ticks(newticks)

    xt,yt=m(-9.10,33.1)
    plt.text(xt,yt,r'($^{\circ}$C)',fontsize=18,\
                    fontname='Times New Roman')
    plt.title(figtitle,fontsize=24,\
                    fontname='Times New Roman')
    
##    ### Export figure and display it
    plt.savefig(figdir+figname+figtype1, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
##    plt.savefig(figdir+figname+figtype2, dpi=300, facecolor='w', edgecolor='w',
##                 transparent=False, bbox_inches='tight', pad_inches=0.1)

    #plt.show()
    plt.close()

