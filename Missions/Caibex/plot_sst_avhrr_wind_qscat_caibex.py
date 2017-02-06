#!/usr/bin/env python

# plot_sst_avhrr_wind_qscat_caibex.py
#
# Plot both wind and SST
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
basemap_resolution = 'l'

# File and directory names
sstdir='/media/ctroupin/Iomega_HDD/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
sstfile='20090827-NAR18_SST-EUR-L2P-sst1nar_noaa18_20090827_desc-v01.nc'
winddir='/home/ctroupin/DataOceano/Satellite/QuikSCAT/netcdf/L2B12v3/data/239/'
windfile='qs_l2b_53067_v3_200908271757.nc'

figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/'
coordfile=sstdir+'coord_NAR_SST.nc'

# Figure extension
figtype1='.png'
figtype2='.eps'

# Colormap
cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps')
cmap1 = colors.LinearSegmentedColormap('my_colormap',cdict,256)
cmap2 = plt.cm.hot

#------------------------------------------------------------------------------------
# Compute min and max values

# 1. SST
vmin1=17.0
vmax1=25.0
norm1 = colors.Normalize(vmin=vmin1,vmax=vmax1)
nlevels1=150
levels2plot1=np.arange(vmin1,vmax1+0.0001,(vmax1-vmin1)/(nlevels1-1))
newticks1=np.arange(np.floor(vmin1),np.ceil(vmax1)+0.0001,1.0)

# 2. Wind 
vmin2 = 0.
vmax2 = 15.
dvar2=2.
bounds2 = np.arange(vmin2,vmax2+0.01,dvar2)
norm2 = colors.Normalize(vmin=vmin2,vmax=vmax2+0.01)
levels2plot2 = np.arange(vmin2,vmax2+0.001,dvar2/10)

#------------------------------------------------------------------------------------
# Region of interest        
lonmin=-13.0
lonmax=-9.0
latmin=29.0
latmax=33.0

# Limits in the matrices (determined emperically) for SST
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
    


# Load coordinates from coordinate file
nc=netcdf.Dataset(coordfile)
lon1 = nc.variables['lon'][imin:imax,jmin:jmax]
lat1 = nc.variables['lat'][imin:imax,jmin:jmax]
nc.close()


#------------------------------------------------------------------------------------
# Load SST
nc=netcdf.Dataset(sstdir+sstfile)
sst = nc.variables['sst'][0,imin:imax,jmin:jmax]
nc.close()
sst = sst.squeeze()

# Mask land values
sst = np.ma.masked_array(sst,sst==sst.fill_value)


#------------------------------------------------------------------------------------
# Load wind
nc=netcdf.Dataset(winddir+windfile)
lon2 = nc.variables['lon'][:]-360
lat2 = nc.variables['lat'][:]
windspeed = nc.variables['retrieved_wind_speed'][:]
winddirection = nc.variables['retrieved_wind_direction'][:]
windtime = nc.variables['time'][:]
nc.close()

lon2=lon2.flatten()
lat2=lat2.flatten()
windspeed=windspeed.flatten()
winddirection=winddirection.flatten()


#------------------------------------------------------------------------------------
# Select sub-region
goodlon = np.nonzero(np.logical_and(lon2<=lonmax+dlon,lon2>=lonmin-dlon))
goodlon = goodlon[0]
lat2 = lat2[goodlon]
lon2 = lon2[goodlon]
windspeed = windspeed[goodlon]
winddirection = -winddirection[goodlon]+90
goodlat = np.nonzero(np.logical_and(lat2<=latmax+dlat,lat2>=latmin-dlat))
goodlat = goodlat[0]
lat2 = lat2[goodlat]
lon2 = lon2[goodlat]
windspeed = windspeed[goodlat]
winddirection = winddirection[goodlat]

#------------------------------------------------------------------------------------
uwind=windspeed*np.cos(np.deg2rad(winddirection))
vwind=windspeed*np.sin(np.deg2rad(winddirection))
uwind=np.ma.masked_where(uwind==uwind.fill_value,uwind)
vwind=np.ma.masked_where(vwind==vwind.fill_value,vwind)
uwind.data[uwind.data==uwind.data.min()]=0
vwind.data[vwind.data==vwind.data.min()]=0
windnorm=np.sqrt(uwind*uwind+vwind*vwind)





# Make the plot
fig=plt.figure()
ax = fig.add_subplot(111)


m.ax=ax


x1,y1 = m(lon1,lat1)
x2,y2 = m(lon2,lat2)
pcm=m.pcolormesh(x1,y1,sst,cmap=cmap1,norm=norm1)

Q=m.quiver(x2,y2,uwind,vwind,windnorm,\
                               units='width',scale=200,width=0.003,norm=norm2,cmap=cmap1)
k = plt.quiverkey(Q, .9, 1.075, 10, r'$10\, ms^{-1}$',labelpos='S',
                                   fontproperties={'weight':'bold','size':'16'},color='k')


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
cbar=fig.colorbar(pcm,cmap=cmap1,orientation='vertical',fraction=0.1,pad=0.02)
newticks1=newticks1.round(2)
cbar.set_ticks(newticks1)

xt,yt=m(-9.10,33.1)
plt.text(xt,yt,r'($^{\circ}$C)',fontsize=18,\
                fontname='Times New Roman')
##plt.title(figtitle,fontsize=24,\
##                fontname='Times New Roman')

##    ### Export figure and display it
##plt.savefig(figdir+figname+figtype1, dpi=300, facecolor='w', edgecolor='w',
 #            transparent=False, bbox_inches='tight', pad_inches=0.1)
##    plt.savefig(figdir+figname+figtype2, dpi=300, facecolor='w', edgecolor='w',
##                 transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
plt.close()

