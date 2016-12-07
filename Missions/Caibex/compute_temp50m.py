#--------------------------------------------------------------------------
# compute_temp50m.py
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

seasoardir = '/home/ctroupin/ULPGC/CAIBEX_campaign/data/Seasor/'
seasoarfile = 'Seasaor2plot.dat'
figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/'
figname='temp50m'
outputdir='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/data/isothermdepth/'
outputfile='temp50m.dat'

# Colormap
cmap=plt.cm.spectral_r

# Load first columns of file
lon,lat,depth,temp=np.loadtxt(seasoardir+seasoarfile,usecols=(0,1,2,3),unpack=True)

##temp=temp[:2000]
##depth=depth[:2000]
depthdiff=depth-50.0
depthsign=depthdiff/abs(depthdiff)
depthzero=np.gradient(depthsign)



goodindices=np.where(abs(depthzero)==1.0)[0]
npoints=len(goodindices)
# Allocate
temp50=np.zeros(npoints)
lon50=np.zeros(npoints)
lat50=np.zeros(npoints)

for t,ind in enumerate(goodindices):

    # Select points for interpolation
    temp2interp=temp[ind-1:ind+2]
    depth2interp=depth[ind-1:ind+2]
    lat2interp=lat[ind-1:ind+2]
    lon2interp=lon[ind-1:ind+2]
    
    # Need to have x increasing before apply interp1d
    index_sort=np.argsort(depth2interp)
    
    # Create interpolator and interpolate
    depth_interpolator = interpolate.interp1d(depth2interp[index_sort],temp2interp[index_sort])
    lon_interpolator = interpolate.interp1d(depth2interp[index_sort],lon2interp[index_sort])
    lat_interpolator = interpolate.interp1d(depth2interp[index_sort],lat2interp[index_sort])
    
    temp50[t]=depth_interpolator(50.0)
    lat50[t]=lat_interpolator(50.0)
    lon50[t]=lon_interpolator(50.0)

##    # plot to check
##    plt.plot(depth2interp,lat2interp,'ko-')
##    plt.plot(50,lat50[t],'ro',color='r',ms=5)
##    plt.show()
##    plt.close()


if (doplot==1):
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
    x,y=m2(lon50,lat50)
    scat=m2.scatter(x,y,c=temp50,edgecolor='none',cmap=cmap)
    m2.drawcoastlines(zorder=4)
    m2.fillcontinents(color='.85',zorder=3)
    cbar2=plt.colorbar(scat,orientation='horizontal',pad=0.02,shrink=0.95,aspect=20)
    cbar2.set_label(r'$(^{\circ}C)$')
    plt.savefig(figdir+figname, facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight')

    #plt.show()
    plt.close()

if dowrite == 1:
    things2write=np.vstack((lon50,lat50,temp50))
    np.savetxt(outputdir+outputfile,things2write.T,delimiter='\t',fmt='%f %f %f')
