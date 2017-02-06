#--------------------------------------------------------------------------
# plot_interpolated_depth_isotherm.py

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


doplot,dowrite=0,1

resdir='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/results/isothermdepth/'
resfile='results.nc'
figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/'
figname='depth_17isotherm_map'
outputdir='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/data/isothermdepth/'
outputfile='isotherm17depth.dat'

fmin=30.
fmax=150.
norm=colors.Normalize(fmin,fmax)
levels2plot=np.arange(fmin,fmax,2)
# Colormap
cmap=plt.cm.spectral_r

# read information from topography file
nc=netcdf.Dataset(resdir+resfile)
lon = nc.variables['x'][:] 
lat = nc.variables['y'][:] 
field = nc.variables['analyzed_field'][:] 
nc.close()

# subregion
lonmin2 = -12.25 
lonmax2 = -9.25 
latmin2 = 29.99 
latmax2 = 31.75 


m2 = Basemap(projection='merc',llcrnrlon=lonmin2,llcrnrlat=latmin2,\
                urcrnrlon=lonmax2,urcrnrlat=latmax2,  \
                lat_ts=0.5*(latmin2+latmax2),\
                resolution='h')
meridians2=np.arange(-12.,-8.5,1.)
parallels2=np.arange(29.5,32.,0.5)
m2.drawmeridians(meridians2,labels=[0, 0, 1, 0],linewidth=0.5,fontsize=18,zorder=2)
m2.drawparallels(parallels2,labels=[0, 1, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

llon,llat=np.meshgrid(lon,lat)
x,y=m2(llon,llat)
x17,y17=m2(lon17,lat17)
#scat=m2.scatter(x17,y17,c=depth17,edgecolor='none',cmap=cmap,norm=norm,zorder=3)
plt.plot(x17,y17,'ko',ms=1)
cont=m2.contourf(x,y,field,levels2plot,cmap=cmap,norm=norm,zorder=2)

m2.drawcoastlines(zorder=4)
m2.fillcontinents(color='.85',zorder=3)

cbar2=plt.colorbar(scat,orientation='horizontal',pad=0.02,shrink=0.95,aspect=20)
cbar2.set_label('(m)')
plt.savefig(figdir+figname, facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight')

#plt.show()
plt.close()


