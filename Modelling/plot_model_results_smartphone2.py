#!/usr/bin/env python
#
# plot_bathymetry_westmedsea.py
#
#
#
# colormap taken from
# http://soliton.vm.bytemark.co.uk/pub/cpt-city/esri/hypsometry/bath/tn/bath_112.png.index.html
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

t0 = time.time()

basemap_resolution = 'l'

#resdir='http://thredds.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmopv2/'
#resfile='latest.nc'

resdir='http://thredds.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmopv2/2014/06/'
resfile='roms_wmopv2_20140622.nc'

figdir='/home/ctroupin/SOCIB/Facilities/Modelling/figures/V14/'
temp_figname='sea_water_potential_temperature_wmopv2'
vel_figname='sea_water_velocity_wmopv2'
salt_figname='sea_water_salinity_wmopv2'
media_type='mobile'

figtype='.png'

# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=-1.0,6.0,37.5,41.5,0.,.0

# Sub-sampling for fields and vectors
MM=1
NN=10           # Only one vector out of 10 is plotted

# temperature
temp_vmin,temp_vmax,temp_dvar = 10.,28.,0.5
# velocity
vel_vmin,vel_vmax,vel_dvar = 0.0,1.0,0.1
# salinity
salt_vmin,salt_vmax,salt_dvar = 36.5,38.5,0.125

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
temp_cmap=plt.cm.RdYlBu_r
salt_cmap=plt.cm.RdYlBu_r
vel_cmap=plt.cm.hot_r
vel_cmap=plt.cm.RdYlBu_r

# Norm and levels to plot
temp_norm = colors.Normalize(vmin=temp_vmin,vmax=temp_vmax+0.001)
temp_levels2plot=np.arange(temp_vmin-0.000,temp_vmax+0.000,temp_dvar)
vel_norm = colors.Normalize(vmin=vel_vmin,vmax=vel_vmax+0.001)
vel_levels2plot=np.arange(vel_vmin-0.000,vel_vmax+0.000,vel_dvar)
salt_norm = colors.Normalize(vmin=salt_vmin,vmax=salt_vmax+0.001)
salt_levels2plot=np.arange(salt_vmin-0.000,salt_vmax+0.000,salt_dvar)

# --------------------------------
# Load results
# --------------------------------

with netcdf.Dataset(resdir+resfile) as nc:

    # Load coordinates
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    lon_uv = nc.variables['lon_uv'][:]
    lat_uv = nc.variables['lat_uv'][:]

    # Select sub-region for temp, salt
    goodlon = np.nonzero(np.logical_and(lon_rho<=lonmax+dlon,lon_rho>=lonmin-dlon))[0]
    goodlat = np.nonzero(np.logical_and(lat_rho<=latmax+dlat,lat_rho>=latmin-dlat))[0]

    lat_rho = lat_rho[goodlat]
    lon_rho = lon_rho[goodlon]
    temp = nc.variables['temp'][:ntimes,goodlat,goodlon]
    salt = nc.variables['salt'][:ntimes,goodlat,goodlon]

    # Select sub-region for u,v 

    goodlon = np.nonzero(np.logical_and(lon_uv<=lonmax+dlon,lon_uv>=lonmin-dlon))[0]
    goodlat = np.nonzero(np.logical_and(lat_uv<=latmax+dlat,lat_uv>=latmin-dlat))[0]
    lat_uv = lat_uv[goodlat]
    lon_uv = lon_uv[goodlon]

    u = nc.variables['u'][:ntimes,goodlat,goodlon]
    v = nc.variables['v'][:ntimes,goodlat,goodlon]
    ocean_time = nc.variables['ocean_time'][:ntimes]
 


# Projection and coordinate grid
llon, llat  = np.meshgrid(lon_rho,lat_rho)
x_rho,y_rho = m(llon, llat)
llon, llat  = np.meshgrid(lon_uv,lat_uv)
x_uv,y_uv = m(llon, llat)

# Prepare time
timeinit = calendar.timegm((date.today().year,date.today().month,date.today().day,0,0,0))
time_model = timeinit+ocean_time

        
# Loop on time
for itime in range(ntimes):


    temp2plot=np.squeeze(temp[itime,:,:])
    salt2plot=np.squeeze(salt[itime,:,:])
    u2plot=np.squeeze(u[itime,:,:])
    v2plot=np.squeeze(v[itime,:,:])
    norm2plot=np.sqrt(u2plot*u2plot+v2plot*v2plot)

    figtime=time.strftime("%Y%m%d_%H%M",time.gmtime(time_model[itime]))
    
    # --------------
    # Make the plots
    #---------------

    # Temperature
    # -------------
    
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.)
    ax = fig.add_subplot(111)
    m.ax=ax
    
    contour=m.contourf(x_rho[::MM],y_rho[::MM],temp2plot[::MM],temp_levels2plot,cmap=temp_cmap,norm=temp_norm,extend='both')

    plt.savefig(figdir+temp_figname+'_'+media_type+'_'+figtime+figtype, dpi=resolution, facecolor='w', edgecolor='w',
      transparent=True, bbox_inches='tight', pad_inches=0.)

    #plt.show()
    plt.close()


    # Velocity
    # -------------
    
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.)
    ax = fig.add_subplot(111)
    m.ax=ax
    
    contour=m.contourf(x_uv[::MM],y_uv[::MM],norm2plot[::MM],vel_levels2plot,cmap=vel_cmap,norm=vel_norm,extend='max')
    Q=plt.quiver(x_uv[::NN, ::NN],y_uv[::NN, ::NN],u2plot[::NN, ::NN],v2plot[::NN, ::NN],\
                    units='width',scale=10,width=0.003,color='k')

    ##plt.colorbar(contour)
    ##
    ##    k = plt.quiverkey(Q, 0.22, 0.9, 0.05, r'$0.05\, ms^{-1}$',labelpos='W',
    ##               fontproperties={'weight':'bold','size':'16'},color='k')

    plt.savefig(figdir+vel_figname+'_'+media_type+'_'+figtime+figtype, dpi=resolution, facecolor='w', edgecolor='w',
      transparent=True, bbox_inches='tight', pad_inches=0.)

    #plt.show()
    plt.close()


    # Salinity
    # -------------
    
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.)
    ax = fig.add_subplot(111)
    m.ax=ax


    contour=m.contourf(x_rho[::MM],y_rho[::MM],salt2plot[::MM],salt_levels2plot,cmap=salt_cmap,norm=salt_norm,extend='both')

    plt.savefig(figdir+salt_figname+'_'+media_type+'_'+figtime+figtype, dpi=resolution, facecolor='w', edgecolor='w',
      transparent=True, bbox_inches='tight', pad_inches=0.)

    #plt.show()
    plt.close()

t1 = time.time()
print str(t1-t0) + ' seconds'
