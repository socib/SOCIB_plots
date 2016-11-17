#!/usr/bin/env python
#
# plot_radar_velocity_vorticity.py
#
# Plot HF-radar velocity with a color corresponding to the
# velocity and maybe the vorticity
#
# ctroupin, July 2014
#----------------------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import time, datetime, calendar

coastdir='/home/ctroupin/DataOceano/Coastlines/'                    # Coastline
coastfile='Ibiza_coastline.dat'
radardir='http://thredds.socib.es/thredds/dodsC/hf_radar/hf_radar_ibiza-scb_codarssproc001/L1/2015/'
radarbasename='dep0001_hf-radar-ibiza_scb-codarssproc001_L1_2015-'
figdir='/home/ctroupin/SOCIB/Facilities/Radar/figures/velocity/2015/'
figbasename='HFradar_vorticity_'
figbasename2='HFradar_velocity_'

doplotvort,doplotvel=1,0

if not(os.path.exists(figdir)):
    os.makedirs(figdir)

figtype='.png'

cmap=plt.cm.RdYlBu_r
basemap_resolution = 'l'
deg2km = 111125.1  

lonmin,lonmax,latmin,latmax,dlon,dlat = 0.5,1.500001,38.35,39.15001,0.2,0.2

lon2plot=np.arange(lonmin,lonmax,dlon)
lat2plot=np.arange(38.4,latmax,dlat)


vmin = -0.5
vmax = 0.5

bounds2 = np.arange(vmin,vmax+0.00001,0.1)
bounds2[len(bounds2)/2]=0
norm2 = colors.Normalize(vmin=vmin,vmax=vmax+0.001)
levels2plot=bounds2



m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            lat_ts=0.5*(latmin+latmax),\
            resolution=basemap_resolution)

i=1
file2load=radardir+radarbasename+str(i).zfill(2)+'.nc'
nc=netcdf.Dataset(file2load)
lon = nc.variables['LON'][:]
lat = nc.variables['LAT'][:]
time0 = nc.variables['time'][:]
u = nc.variables['U'][:]
v = nc.variables['V'][:]
#norm = nc.variables['WSPE'][:]
nc.close()

# Compute vorticity
lon,lat=np.meshgrid(lon,lat)
dlonx,dlony=np.gradient(lon)
dlatx,dlaty=np.gradient(lat)
dx=dlony*np.cos(np.deg2rad(lat))*deg2km
dy=dlatx*deg2km
f=4*np.pi*np.sin(np.deg2rad(lat))/86400

NN = 1
x,y = m(lon, lat)

# ----------------------------------------------------------------------------
# Load coastline extracted from
# http://www.ngdc.noaa.gov/mgg_coastline/

valexc=-999.
lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
lonc=np.ma.masked_where(lonc==valexc,lonc)
latc=np.ma.masked_where(latc==valexc,latc)
xc,yc=m(lonc,latc)

count=0
for i in range(3,4):

    file2load=radardir+radarbasename+str(i).zfill(2)+'.nc'
    print 'Working on file '+os.path.basename(file2load)

    nc=netcdf.Dataset(file2load)
    time0 = nc.variables['time'][:]
    u = nc.variables['U'][:]
    v = nc.variables['V'][:]
    nc.close()

    # Mask 
    u = np.ma.masked_array(u,np.isnan(u))
    v = np.ma.masked_array(v,np.isnan(v))
    
    for t,k in enumerate(time0):
        
        
        u2plot = u[t,:,:]
        v2plot = v[t,:,:]
        time2plot= time0[t]
        norm2plot= np.sqrt(u2plot*u2plot*v2plot*v2plot)

        # Vorticity
        dux,duy=np.gradient(np.squeeze(u2plot))
        dvx,dvy=np.gradient(np.squeeze(v2plot))
        vort=(dvy/dx-dux/dy)/f

        
        figsuffix = time.strftime('%Y_%m_%d_%H_%M_%S', time.gmtime(time2plot))
        figtitle = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time2plot))


        if doplotvort==1:
            # Make the plot
            fig=plt.figure()
            ax = fig.add_subplot(111)
            
            m.ax=ax
                
            Q=plt.quiver(x[::NN, ::NN],y[::NN, ::NN],u2plot[::NN, ::NN],v2plot[::NN, ::NN],vort[::NN, ::NN],\
                         units='width',scale=7.5,width=0.003,color='k',zorder=2,cmap=cmap,norm=norm2)

            #pcm=m.contourf(x,y,vort,levels2plot,zorder=1,cmap=cmap,norm=norm2,extend='both')
            k = plt.quiverkey(Q, 0.9, 0.95, 0.25, r'$0.25\, ms^{-1}$',labelpos='W',
                       fontproperties={'weight':'bold','size':'16'},color='k')

            m.drawparallels(lat2plot, linewidth=0.5,labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=1)
            m.drawmeridians(lon2plot, linewidth=0.5,labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16,zorder=1)


            m.plot(xc,yc,lw=0.5,color='k')
            m.fillcontinents(color='w', ax=ax,zorder=3)

            # Add the colorbar
            cbar=fig.colorbar(Q,cmap=cmap,norm=norm2,orientation='vertical',pad=0.025,aspect=15,extend='both')
            cbar.set_label(r'$\xi/f$',fontname='Times New Roman',fontsize=18,rotation=0)
            cbar.set_ticks(bounds2)
            cbar.set_clim((vmin,vmax,))
            cbar.ax.set_yticklabels(bounds2,fontname='Times New Roman',fontsize=16)
            cbar.set_cmap(plt.cm.spectral)
            cbar.set_norm(norm2)
            plt.title(figtitle,fontname='Times New Roman',fontsize=24,color='k')

            count+=1
            figname=figdir+figbasename+str(count).zfill(3)+figtype
            plt.savefig(figname, dpi=150, facecolor='w', edgecolor='w',
                    transparent=False, bbox_inches='tight', pad_inches=0.1)
            
            #plt.show()
            plt.close()

        if doplotvel==1:
            # Make the plot
            fig=plt.figure()
            ax = fig.add_subplot(111)
            
            m.ax=ax
                
            Q=plt.quiver(x[::NN, ::NN],y[::NN, ::NN],u2plot[::NN, ::NN],v2plot[::NN, ::NN],norm2plot,\
                         units='width',scale=7.5,width=0.003,color='k',zorder=2,cmap=cmap)

            
            k = plt.quiverkey(Q, 0.9, 0.95, 0.25, r'$0.25\, ms^{-1}$',labelpos='W',
                       fontproperties={'weight':'bold','size':'16'},color='k')

            m.drawparallels(lat2plot, linewidth=0.5,labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=2)
            m.drawmeridians(lon2plot, linewidth=0.5,labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16,zorder=2)


            m.plot(xc,yc,lw=0.5,color='k')
            m.fillcontinents(color='w', ax=ax,zorder=3)

            plt.title(figtitle,fontname='Times New Roman',fontsize=24,color='k')

            count+=1
            figname=figdir+figbasename2+figsuffix+figtype
##            plt.savefig(figname, dpi=75, facecolor='w', edgecolor='w',
##                    transparent=False, bbox_inches='tight', pad_inches=0.1)
            
            plt.show()
            plt.close()
