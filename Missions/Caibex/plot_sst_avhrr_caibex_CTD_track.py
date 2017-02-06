#!/usr/bin/env python

# plot_sst_avhrr_caibex_CTD_track.py
#
# Plot the SST 
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
CTDdir='/home/ctroupin/ULPGC/CAIBEX_campaign/data/CTD/'
CTDfile='CTD_coordinates_transect.dat'
datadir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/4paperCAIBEX/'
databasename='2009*.nc'
figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'

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
newticks=newticks.round(2)
    
# Region of interest        
lonmin=-13.0
lonmax=-8.9
latmin=29.0
latmax=33.0

# Limits in the matrices (determined emperically)
imin,imax,jmin,jmax=450,800,1500,1800

# Spacing between x/y labels
dlon = 1
dlat = 1

# Unit of the variable to plot
unitname = ' '


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

# Load CTD coordinates
coordCTD=np.loadtxt(CTDdir+CTDfile)
latCTD=coordCTD[:,0]+coordCTD[:,1]/60.+coordCTD[:,2]/3600.
lonCTD=-1*(coordCTD[:,3]+coordCTD[:,4]/60.+coordCTD[:,5]/3600.)
    
# Load coordinates from coordinate file
nc=netcdf.Dataset(coorddir+coordfile)
lon = nc.variables['lon'][imin:imax,jmin:jmax]
lat = nc.variables['lat'][imin:imax,jmin:jmax]
nc.close()

##print lon.min()
##print lon.max()
##print lat.min()
##print lat.max()


meridians=np.arange(lonmin,lonmax+0.001,dlon)
parallels=np.arange(latmin,latmax+0.001,dlat)

# Projection
x,y = m(lon,lat)
xCTD,yCTD = m(lonCTD,latCTD)


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

    # Mask land values
    field = np.ma.masked_array(field,field==field.fill_value)
    

    # Make the plot
    fig=plt.figure()
    ax = fig.add_subplot(111)
    m.ax=ax
   
    pcm=m.pcolormesh(x,y,field,cmap=cmap,norm=norm)

    # Add CTD coordinates
    m.plot(xCTD,yCTD,'ko',ms=2)

    # Add grid, coastline and continent
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='.15', ax=ax,zorder=3)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

    # Add the colorbar
    cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
    cbar.set_ticks(newticks)

    xt,yt=m(-9.10,33.1)
    plt.text(xt,yt,r'($^{\circ}$C)',fontsize=18)
    plt.title(figtitle,fontsize=24)
    
##    ### Export figure and display it
    plt.savefig(figdir+figname+figtype1, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
##    plt.savefig(figdir+figname+figtype2, dpi=300, facecolor='w', edgecolor='w',
##                 transparent=False, bbox_inches='tight', pad_inches=0.1)
##
##    plt.show()
    plt.close()

