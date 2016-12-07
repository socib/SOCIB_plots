#!/usr/bin/env python

# CapeGhir_plot_mean_ADT.py
#
# Compute the average over the period of the cruise
# plot the ADT and the velocity around Cape Ghir
#
# -------------------------------------------------------------


import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
origin = 'lower'

# Plot the gridded field obtained from AVISO FTP

fieldtype='SLA'
os.system('clear')
basemap_resolution = 'i'
filedir='/home/ctroupin/DataOceano/AVISO/Canary/data/'+fieldtype+'/'
figdir='/home/ctroupin/DataOceano/AVISO/Canary/figures/'+fieldtype+'/'
figtype='.png'

cmap=plt.cm.RdBu_r


lonmin=-13
lonmax=-9
latmin=29
latmax=33
dlon = 1
dlat = 1

vmin = -7.50
vmax = 7.5
bounds2 = np.arange(vmin,vmax+0.01,5)
norm = colors.Normalize(vmin=vmin,vmax=vmax+0.01)

nlevels=50
dvar=(np.real(vmax)-np.real(vmin))/np.real(nlevels)
levels2plot=np.arange(vmin-0.01,vmax+0.01,dvar)
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            resolution=basemap_resolution)

NN=1

filelist=sorted(glob.glob(filedir + 'dt_upd_global_merged_msla_h_*.nc') )

# Load info in the first file
firstfile=filelist[0]
nc=netcdf.Dataset(firstfile)
lon = nc.variables['NbLongitudes'][:]-360.0
lat = nc.variables['NbLatitudes'][:]
nc.close()

# Select sub-region (faster for plotting)
goodlon = np.nonzero(np.logical_and(lon<=lonmax+2*dlon,lon>=lonmin-2*dlon))
goodlat = np.nonzero(np.logical_and(lat<=latmax+2*dlat,lat>=latmin-2*dlat))
goodlat = goodlat[0][:]
goodlon = goodlon[0][:]

lat = lat[goodlat]
lon = lon[goodlon]

llon, llat  = np.meshgrid(lon,lat)

ntimes=len(filelist)
nlat=len(lat)
nlon=len(lon)

# Allocate
field3D=np.zeros((ntimes,nlat,nlon,))
u3D=np.zeros((ntimes,nlat,nlon,))
v3D=np.zeros((ntimes,nlat,nlon,))

figname=figdir+'meanSLA'+figtype

i=0
for infile in filelist:
    filename= os.path.basename(infile)
    print('Working on file ' + filename)
    

    # Now prepare velocity file    
    filedate=filename.rsplit("_")[6]
    velocityfile=glob.glob(filedir + 'dt_upd_global_merged_msla_uv_'+filedate+'*')


    # Load ADT
    nc=netcdf.Dataset(filedir+filename)
    field = nc.variables['Grid_0001'][:]
    field = field.T         # Transpose
    nc.close()

    field = field[goodlat]
    field = field[:,goodlon]

    #Mask the field
    valex = -99;
    field = np.ma.masked_where(field>1000, field)

    # Write in 3D matrix
    field3D[i,:,:]=field
    

    # Load velocity components
    nc=netcdf.Dataset(velocityfile[0])
    lonUV = nc.variables['NbLongitudes'][:]-360
    latUV = nc.variables['NbLatitudes'][:]
    u = 0.01*nc.variables['Grid_0001'][:].T
    v = 0.01*nc.variables['Grid_0002'][:].T
    nc.close()

    u = u[goodlat]
    u = u[:,goodlon]
    v = v[goodlat]
    v = v[:,goodlon]

    u3D[i,:,:]=u
    v3D[i,:,:]=v

    i+=1

# Compute time-averaged fields
meanfield=np.mean(field3D,axis=0)
meanu=np.mean(u3D,axis=0)
meanv=np.mean(v3D,axis=0)

meanfield = np.ma.masked_where(meanfield>100, meanfield)
meanu = np.ma.masked_where(meanu>100, meanu)
meanv = np.ma.masked_where(meanv>100, meanv)

# Make the plot

fig=plt.figure()
ax = fig.add_subplot(111)
m.ax=ax
x,y = m(llon, llat)
contour=m.contourf(x,y,meanfield,levels2plot,cmap=cmap,norm=norm,extend='both')
#contour2=m.contour(x,y,field,levels2plot,linewidths=0.5,colors='k',zorder=1,origin=origin,)

Q=plt.quiver(x[::NN, ::NN],y[::NN, ::NN],meanu[::NN, ::NN],meanv[::NN, ::NN],\
            units='width',scale=2.5,width=0.003,color='k')
k = plt.quiverkey(Q, 0.5, 0.05, 0.1, r'$0.1\, ms^{-1}$',labelpos='W',
           fontproperties={'weight':'bold','size':'16'},color='k',zorder=3)  

m.drawcoastlines(ax=ax)
m.fillcontinents(zorder=2)

m.drawparallels(np.arange(latmin,latmax+0.0001,dlat), linewidth=0.5,
                    labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16)
m.drawmeridians(np.arange(lonmin,lonmax+0.0001,dlon), linewidth=0.5,
                    labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16)
##    
# Add the colorbar
cbar=fig.colorbar(contour,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
cbar.set_label('(cm)',fontname='Times New Roman',fontsize=18)
cbar.set_ticks(bounds2)
cbar.ax.set_yticklabels(bounds2,fontname='Times New Roman',fontsize=16)

#plt.savefig(figname, dpi=300, facecolor='w', edgecolor='w',
#         transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()

