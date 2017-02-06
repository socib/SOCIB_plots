#--------------------------------------------------------------------------
# plot_domain_caibex2014.m
#
#
# ctroupin, July 2012
# Adapted September 2012
# Python version January 2014
#--------------------------------------------------------------------------

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
from matplotlib.path import Path
import matplotlib.patches as patches

# files and directories

figdir = '/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/' 
figname = '1_CaibexBathymetry30sec_V2' 
griddir = '/home/ctroupin/DataOceano/Topo_GEBCO/'
gridfile = 'gebco_08_-30_20_0_40.nc' 
CTDdir='/home/ctroupin/ULPGC/CAIBEX_campaign/data/CTD/'
CTDfile='CTD_coordinates_transect.dat'
seasoardir = '/home/ctroupin/ULPGC/CAIBEX_campaign/data/Seasor/'
seasoarfile = 'Seasaor2plot.dat'

# cape ghir coordinates
lat_g = 30.+(37.+49./60.)/60. 
lon_g = -(9.+(53.+20./60.)/60.) 

contourlevels = np.array([-500., -1000., -2000., -3000., -4000.])

meridians1=np.arange(-18.,-5.,3.)
parallels1=np.arange(28.,37.,2.)
meridians2=np.arange(-12.,-8.5,1.)
parallels2=np.arange(29.5,32.,0.5)

cmap=plt.cm.gist_earth

# Prepare colormaps
cdict= gmtColormap.gmtColormap('bath_112','/home/ctroupin/Software/Python/GMT_colormaps')
cmap_ocean = colors.LinearSegmentedColormap('my_colormap',cdict,256)
cdict= gmtColormap.gmtColormap('spain','/home/ctroupin/Software/Python/GMT_colormaps')
cmap_land = colors.LinearSegmentedColormap('my_colormap',cdict,256)

levels2plot_ocean = np.arange(-6000,0.1,100)
levels2plot_ocean2 = np.arange(-3500,0.1,100)
levels2plot_land = np.arange(0,3800,100)

norm1 = colors.Normalize(vmin=-6000,vmax=0)
norm2 = colors.Normalize(vmin=-3500,vmax=0)
#--------------------------------------------------------------------------

if not(os.path.exists(figdir)):
    os.makedirs(figdir)

# read information from topography file
nc=netcdf.Dataset(griddir+gridfile)
lontopomin,lontopomax = nc.variables['x_range'][0], nc.variables['x_range'][1] 
lattopomin,lattopomax = nc.variables['y_range'][0], nc.variables['y_range'][1]
dlon = nc.variables['spacing'][0] 
dlat = nc.variables['spacing'][1] 
bathy = nc.variables['z'][:] 
nc.close()

# Create the vectors
lontopo = np.arange(lontopomin,lontopomax,dlon)
lattopo = np.arange(lattopomin,lattopomax,dlat)
bathy = np.reshape(bathy,(len(lattopo),len(lontopo)))
bathy = np.flipud(bathy)

##### Graphical check
##plt.pcolormesh(lontopo,lattopo,bathy)
##plt.show()
##plt.close()


# Sub-sampliong (if necessary)
NN = 5 
bathy = bathy[0::NN,0::NN]
lontopo = lontopo[0::NN] 
lattopo = lattopo[0::NN]

# Load CTD coordinates
coordCTD=np.loadtxt(CTDdir+CTDfile)
latCTD=coordCTD[:,0]+coordCTD[:,1]/60.+coordCTD[:,2]/3600.
lonCTD=-1*(coordCTD[:,3]+coordCTD[:,4]/60.+coordCTD[:,5]/3600.)

# Load SeaSoar tracks
#-----------------------

lonseasoar,latseasoar = np.loadtxt(seasoardir+seasoarfile,usecols=(0,1),unpack='True')


# prepare projection
#-------------------

# domain of interest
lonmin = -19 
lonmax = -5. 
latmin = 27. 
latmax = 36.

# subregion
lonmin2 = -12.25 
lonmax2 = -9.25 
latmin2 = 29.99 
latmax2 = 31.75 

m1 = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                lat_ts=0.5*(latmin+latmax),\
                resolution='h')

m2 = Basemap(projection='merc',llcrnrlon=lonmin2,llcrnrlat=latmin2,\
                urcrnrlon=lonmax2,urcrnrlat=latmax2,  \
                lat_ts=0.5*(latmin2+latmax2),\
                resolution='f')

##
### Remove bad point
###-----------------
goodlon = np.nonzero(np.logical_and(lontopo<=lonmax+dlon,lontopo>=lonmin-dlon))
goodlat = np.nonzero(np.logical_and(lattopo<=latmax+dlat,lattopo>=latmin-dlat))
goodlon,goodlat=goodlon[0],goodlat[0]
lontopo = lontopo[goodlon] 
lattopo = lattopo[goodlat]
bathy = bathy[goodlat,:]
bathy = bathy[:,goodlon]

goodlon = np.nonzero(np.logical_and(lontopo<=lonmax2+dlon,lontopo>=lonmin2-dlon))
goodlat = np.nonzero(np.logical_and(lattopo<=latmax2+dlat,lattopo>=latmin2-dlat))
goodlon,goodlat=goodlon[0],goodlat[0]
lontopo2 = lontopo[goodlon] 
lattopo2 = lattopo[goodlat]
bathy2 = bathy[goodlat,:]
bathy2 = bathy2[:,goodlon]


### Graphical check
##plt.pcolormesh(lontopo,lattopo,bathy)
##plt.show()
##plt.close()

# Prepare text and projection
xt1,yt1=m1(-18.5,35.5)
xt2,yt2=m1(-18.5,27.2)
xt3,yt3=m1(-18.5,33.8)
xt4,yt4=m1(-8.,27.2)

xG,yG=m1(-9.5,30.55) 
xS,yS=m1(-9.5,31.4) 
xJ,yJ=m1(-12.9,27.55) 
xB,yB=m1(-8.8,32.7) 
xG2,yG2=m2(-9.8,30.6) 
xS2,yS2=m2(-9.75,31.4) 

# ---------------
# Make the plot
# ---------------

llontopo,llattopo=np.meshgrid(lontopo,lattopo)
xtopo,ytopo=m1(llontopo,llattopo)

fig=plt.figure(num=None, figsize=(14,5.5), facecolor='w', edgecolor='k')


# 1. Large subdmain
ax = fig.add_subplot(121)

contour1=m1.contourf(xtopo,ytopo,bathy,levels2plot_ocean,edgecolor=None,cmap=cmap_ocean,norm=norm1,alpha=0.85)
m1.contourf(xtopo,ytopo,bathy,levels2plot_land,edgecolor=None,cmap=cmap_land,zorder=3)

m1.drawmeridians(meridians1,labels=[0, 0, 1, 0],linewidth=0.5,fontsize=18,zorder=2)
m1.drawparallels(parallels1,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)
m1.drawcoastlines(ax=ax,zorder=4)
#m1.fillcontinents(color='.85', ax=ax,zorder=3)

# Add dots at capes
xcape,ycape=m1(-(9.+(51./60.)),(31.+24./60.))
m1.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)  # Sim
xcape,ycape=m1(-(9.+(16./60.)),(32.+33./60.))
m1.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)  # Bedouzza
xcape,ycape=m1(-(12.+(56./60.)),(27.+56./60.))
m1.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)  # Jubi
xcape,ycape=m1(-(14.+(24./60.)),(26.+12./60.))
m1.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)  # Bojador
xcape,ycape=m1(-(17.+(3./60.)),(20.+58./60.))
m1.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)  # Blanc
xcape,ycape=m1(lon_g,lat_g)
m1.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)

# Add cape names
plt.text(xt1,yt1,'NORTH ATLANTIC OCEAN',color='w')
plt.text(xt2,yt2,'Canary Islands',color='w') 
plt.text(xt3,yt3,'Madeira',color='w') 
plt.text(xG,yG,'C. Ghir',zorder=4) 
plt.text(xS,yS,'C. Sim',zorder=4) 
plt.text(xJ,yJ,'C. Juby',zorder=4) 
plt.text(xB,yB,'C. Bedouzza',zorder=4) 
plt.text(xt4,yt4,'AFRICA')

# Add rectangle around region of interest

lonrect=np.array((lonmin2,lonmax2,lonmax2,lonmin2,lonmin2))
latrect=np.array((latmin2,latmin2,latmax2,latmax2,latmin2))

xrect,yrect=m1(lonrect,latrect)
m1.plot(xrect,yrect,'k-')

# Prepare colorbar
cbar=fig.colorbar(contour1,orientation='horizontal',pad=0.02,aspect=15)
cbar.set_ticks([])
cbar.set_ticks(np.arange(-5000,0,1000))
cbar.set_label('(m)')

# 2. Smaller subdmain
ax = fig.add_subplot(122)

llontopo2,llattopo2=np.meshgrid(lontopo2,lattopo2)
xtopo2,ytopo2=m2(llontopo2,llattopo2)
xCTD,yCTD = m2(lonCTD,latCTD)
xseasoar,yseasoar = m2(lonseasoar,latseasoar)

contourf2=m2.contourf(xtopo2,ytopo2,bathy2,levels2plot_ocean2,edgecolor=None,cmap=cmap_ocean,norm=norm2,alpha=0.85)
m2.contourf(xtopo2,ytopo2,bathy2,levels2plot_land,edgecolor=None,cmap=cmap_land,zorder=3)

m2.drawmeridians(meridians2,labels=[0, 0, 1, 0],linewidth=0.5,fontsize=18,zorder=2)
m2.drawparallels(parallels2,labels=[0, 1, 0, 0],linewidth=0.5,fontsize=18,zorder=2)
m2.drawcoastlines(ax=ax,zorder=4)

# Add dots at capes
xcape,ycape=m2(-(9.+(51./60.)),(31.+24./60.))
m2.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)  # Sim
xcape,ycape=m2(lon_g,lat_g)
m2.plot(xcape,ycape,'ko',ms=5,color='k',zorder=4)

# Add cape names
plt.text(xG2,yG2,'C. Ghir',zorder=4) 
plt.text(xS2,yS2,'C. Sim',zorder=4)

# Add CTD and SeaSoar positions
m2.plot(xCTD,yCTD,'ko',ms=5,label='CTD',zorder=4)
m2.plot(xseasoar,yseasoar,'--',ms=0.5,label='SeaSoar',color='0.2',zorder=3)
plt.legend(loc=4)
cbar2=fig.colorbar(contourf2,orientation='horizontal',pad=0.02,aspect=15)
cbar2.set_ticks([])
cbar2.set_ticks(np.arange(-3000,0,500))
cbar2.set_label('(m)')

plt.text(xCTD[0]+10000,yCTD[0],'T1',ha='center',fontweight='bold')
plt.text(xCTD[-1]+10000,yCTD[-1]-15000,'T13',ha='center',fontweight='bold')

xt,yt=m2(-12.,30.025)
plt.text(xt,yt,1,ha='center',fontweight='bold')
xt,yt=m2(-11.75,31.5)
plt.text(xt,yt,2,ha='center',fontweight='bold')
xt,yt=m2(-11.45,30.4)
plt.text(xt,yt,3,ha='center',fontweight='bold')
xt,yt=m2(-11.19,31.3)
plt.text(xt,yt,4,ha='center',fontweight='bold')
xt,yt=m2(-10.87,30.14)
plt.text(xt,yt,5,ha='center',fontweight='bold')
xt,yt=m2(-10.59,31.5)
plt.text(xt,yt,'6-7',ha='center',fontweight='bold')

plt.savefig(figdir+figname, facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight')
plt.show()
plt.close()
