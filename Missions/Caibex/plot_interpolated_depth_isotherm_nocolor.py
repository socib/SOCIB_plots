#--------------------------------------------------------------------------
# plot_interpolated_depth_isotherm_nocolor
#
#
#
# ctroupin, January 2014
#--------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
from scipy import interpolate
from mpl_toolkits.basemap import Basemap

doplot,dowrite=0,1

resdir1='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/results/isothermdepth/'
resfile1='temp50m.nc'
resdir2='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/results/isothermdepth/'
resfile2='results.nc'
figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/'
figname='depth_17isotherm_map'
outputdir='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/data/isothermdepth2/'
outputfile='isotherm17depth.dat'
datafile1='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/data/isothermdepth/temp50m.dat'

# Colormap
cmap=plt.cm.spectral_r

fmin1=16.0
fmax1=21.
norm1=colors.Normalize(fmin1,fmax1)
levels2plot1=np.arange(fmin1,fmax1,0.25)

fmin2=30.
fmax2=150.
norm2=colors.Normalize(fmin2,fmax2)
levels2plot2=np.arange(fmin2,fmax2,5)



# Define domain and projection
lonmin2 = -12.25 
lonmax2 = -9.25 
latmin2 = 29.99 
latmax2 = 31.75 

m2 = Basemap(projection='merc',llcrnrlon=lonmin2,llcrnrlat=latmin2,\
                urcrnrlon=lonmax2,urcrnrlat=latmax2,  \
                lat_ts=0.5*(latmin2+latmax2),\
                resolution='h')

londata,latdata,data=np.loadtxt(datafile1,usecols=(0,1,2),unpack=True)

# read information from topography file
nc=netcdf.Dataset(resdir1+resfile1)
lon1 = nc.variables['x'][:] 
lat1 = nc.variables['y'][:] 
field1 = nc.variables['analyzed_field'][:]
error1 = nc.variables['error_field'][:]
nc.close()

nc2=netcdf.Dataset(resdir2+resfile2)
lon2 = nc2.variables['x'][:] 
lat2 = nc2.variables['y'][:] 
field2 = nc2.variables['analyzed_field'][:]
error2 = nc2.variables['error_field'][:]
nc2.close()

# Mask field where error higher than 0.3 (or other value)
field1 = np.ma.masked_where(error1>=0.25,field1)
field2 = np.ma.masked_where(error2>=0.25,field2)



llon1,llat1=np.meshgrid(lon1,lat1)
x1,y1=m2(llon1,llat1)
llon2,llat2=np.meshgrid(lon2,lat2)
x2,y2=m2(llon2,llat2)
xdata,ydata=m2(londata,latdata)

# Make the plot

fig=plt.figure(num=None, figsize=(14,5.5), facecolor='w', edgecolor='k')

# Temperature at 50m
# -------------------

ax = fig.add_subplot(121)
cont1=m2.contour(x1,y1,field1,levels2plot1,zorder=2,colors='k')
#scat=m2.scatter(xdata,ydata,c=data,edgecolor='none',zorder=3)
plt.clabel(cont1,levels2plot1[0::2],inline=1,fmt='%1.1f',fontsize=14)
m2.drawcoastlines(zorder=4)
m2.fillcontinents(color='.85',zorder=3)
#cbar2=plt.colorbar(cont1,orientation='horizontal',pad=0.02,shrink=0.95,aspect=20)

meridians2=np.arange(-12.,-8.5,0.5)
parallels2=np.arange(29.5,32.,0.5)
m2.drawmeridians(meridians2,labels=[0, 0, 1, 0],linewidth=0.5,fontsize=18,zorder=2)
m2.drawparallels(parallels2,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

ax = fig.add_subplot(122)
#x17,y17=m2(lon17,lat17)
#scat=m2.scatter(x17,y17,c=depth17,edgecolor='none',cmap=cmap,norm=norm,zorder=3)
#plt.plot(x17,y17,'ko',ms=1)
cont2=m2.contour(x2,y2,field2,levels2plot2,colors='k',zorder=2)
plt.clabel(cont2,levels2plot2[0::2],inline=1,fmt='%1.0f',fontsize=14)
m2.drawcoastlines(zorder=4)
m2.fillcontinents(color='.85',zorder=3)
meridians2=np.arange(-12.,-8.5,0.5)
parallels2=np.arange(29.5,32.,0.5)
m2.drawmeridians(meridians2,labels=[0, 0, 1, 0],linewidth=0.5,fontsize=18,zorder=2)
#m2.drawparallels(parallels2,labels=[0, 1, 0, 0],linewidth=0.5,fontsize=18,zorder=2)



##cbar2=plt.colorbar(cont,orientation='horizontal',pad=0.02,shrink=0.95,aspect=20)
##cbar2.set_label('(m)')
plt.savefig(figdir+figname+'_V2', facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight')

plt.show()
plt.close()


