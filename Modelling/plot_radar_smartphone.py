#!/usr/bin/env python
#
# plot_radar_smartphone.py
#
#
#
# ---------------------------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from mpl_toolkits.basemap import Basemap
from matplotlib import rc
from matplotlib  import colors
import time, calendar
from datetime import date
import sys


basemap_resolution = 'l'

resdir='http://thredds.socib.es/thredds/dodsC/hf_radar/hf_radar_ibiza-scb_codarssproc001/L1/'
resfile='dep0001_hf-radar-ibiza_scb-codarssproc001_L1_latest.nc'
figdir='/home/ctroupin/SOCIB/Facilities/Radar/figures/'
radar_figname='sea_water_velocity_radar'
media_type='mobile'

figtype='.png'

# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=0.575,1.35,38.375,39.05,0.0,0.

# Sub-sampling for fields and vectors
MM=1
NN=1           # Only one vector out of 10 is plotted

# velocity
vel_vmin,vel_vmax,vel_dvar = 0.0,1.0,0.1

resolution=75

ntimes=25        # number of time steps to plot (25 to get one day)

# ----------------------------------------------------------------------------

# Create directory if necessary
if not os.path.exists(figdir):
    os.makedirs(figdir)
        

# Initialise basemap projection
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            lat_ts=0.5*(latmin+latmax),\
            resolution=basemap_resolution)

# Colormap
vel_cmap=plt.cm.RdYlBu_r

# Norm and levels to plot
vel_norm = colors.Normalize(vmin=vel_vmin,vmax=vel_vmax+0.001)
vel_levels2plot=np.arange(vel_vmin-0.000,vel_vmax+0.000,vel_dvar)


# --------------------------------
# Load results
# --------------------------------

with netcdf.Dataset(resdir+resfile) as nc:

    # Load coordinates
    lon = nc.variables['LON'][:]
    lat = nc.variables['LAT'][:]
    
    u = nc.variables['U'][-ntimes:,:,:]
    v = nc.variables['V'][-ntimes:,:,:]
    ocean_time = nc.variables['time'][-ntimes:]
 


# Projection and coordinate grid
llon, llat  = np.meshgrid(lon,lat)
x,y = m(llon, llat)

       
# Loop on time
for itime in range(ntimes):


    u2plot=np.squeeze(u[itime,:,:])
    v2plot=np.squeeze(v[itime,:,:])
    norm2plot=np.sqrt(u2plot*u2plot+v2plot*v2plot)

    figtime=time.strftime("%Y%m%d_%H%M",time.gmtime(ocean_time[itime]))
    
   
    
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.)
    ax = fig.add_subplot(111)
    m.ax=ax
    
    contour=m.contourf(x[::MM],y[::MM],norm2plot[::MM,::MM],vel_levels2plot,cmap=vel_cmap,norm=vel_norm,extend='max')
    Q=plt.quiver(x[::NN, ::NN],y[::NN, ::NN],u2plot[::NN, ::NN],v2plot[::NN, ::NN],\
                    units='width',scale=10,width=0.003,color='k')

    ##plt.colorbar(contour)
    ##
    ##    k = plt.quiverkey(Q, 0.22, 0.9, 0.05, r'$0.05\, ms^{-1}$',labelpos='W',
    ##               fontproperties={'weight':'bold','size':'16'},color='k')

##    plt.savefig(figdir+radar_figname+'_'+media_type+'_'+figtime+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)

    plt.show()
    plt.close()

