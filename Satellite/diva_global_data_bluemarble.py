#!/usr/bin/env python

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
from matplotlib.path import Path
import matplotlib.patches as patches


#------------------------------------------------
# Plot the finite-element mesh and the topography
#------------------------------------------------

os.system('clear')
basemap_resolution = 'l'

datadir='/home/ctroupin/DataOceano/AVISO/netcdf_tracks/global/ascii/'
databasename='nrt_global_*20130801.dat'
resdir='/home/ctroupin/DIVA/Global/results/time_weight_5/'
filebasename='nrt_global_merged_sla_vfec_2013'
meshdir='/home/ctroupin/DIVA/Global/mesh/L3/'
figdir='/home/ctroupin/DIVA/Global/figures/results/time_weight_5/'
figtype='.png'

filelist=sorted(glob.glob(datadir+databasename))


# Create figure directory if neccesary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

# Select colormaps
cmap=plt.cm.spectral

vmin = -20.
vmax = 20.
bounds2 = np.arange(vmin,vmax+0.01,5)
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.01)

nlevels=50
dvar=(np.real(vmax)-np.real(vmin))/np.real(nlevels)
levels2plot=np.arange(vmin-0.01,vmax+0.01,dvar)

# Region of interest
# (could be read from GridInfo.dat file)
lonmin=-180.0
lonmax=180.0
latmin=-90.0
latmax=90.0

dlat=30
dlon=60

lonmap=0
latmap=20

#m = Basemap(resolution=basemap_resolution,projection='ortho',lat_0=latmap,lon_0=latmap)

# Create basemap then used for the plot
m = Basemap(projection='ortho', lat_0 = 15, lon_0 = 20,
                  resolution = 'c')
m2 = Basemap(projection='ortho', lat_0 = 15, lon_0 = 90,
                  resolution = 'c')
#------------------------------------------------------------------------------

fig=plt.figure()
ax = fig.add_subplot(111)
m.ax=ax
    
for infile in filelist:
    filename= os.path.basename(infile)
    print('Working on file ' + filename)
    figname1=figdir+filename[0:-3]+'_atlantic'+figtype

    lon,lat,field=np.loadtxt(infile,usecols=(0,1,2),unpack='True')
    field*=100.
##    # Load results
##    nc=netcdf.Dataset(infile)
##    lon = nc.variables['x'][:]
##    lat = nc.variables['y'][:]
##    field = nc.variables['analyzed_field'][:]
##    nc.close()
##
##    #Mask the field
##    valex = -99;
##    field = np.ma.masked_where(field==valex, field)
##    field = field*100.0

    #----------------    
    # Make the plots
    #----------------

    

    x,y = m(lon, lat)
    scat=m.scatter(x,y,c=field,s=5,edgecolor='none',norm=norm,cmap=cmap,zorder=3)

    m.bluemarble()
    #m.drawcountries()

    
plt.savefig(figname1, dpi=150, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()

