#--------------------------------------------------------------------------
# plot_CTD_transect.py
#
# plot the interpolated results over the CTD tracks
#
# Same kind of subfigure than Pablo
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
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


doplot,dowrite=0,1
valex=-99.0
CTDdir='/home/ctroupin/ULPGC/CAIBEX_campaign/data/CTD/'
CTDfile='CTD_coordinates_transect.dat'
resdir='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/results/transectCTD/L10SNR1icoord80/'
figdir='/home/ctroupin/ULPGC/CAIBEX_campaign/diva/figures/analysisCTDVertical_python/'
resfiles=np.array(('tempL10SNR1icoord80.nc','psalL10SNR1icoord80.nc'))


depthmin=400
depth2plot=range(-depthmin,0,50)
depth2plot2=range(depthmin,0,-50)

levels2plot=np.arange(10,23,0.5)
levels2plot2=np.arange(35,37.8,0.1)


stationlabels=('13','','11','','9','','7','','5','','3','','1')


# Load CTD coordinates
coordCTD=np.loadtxt(CTDdir+CTDfile)
latCTD=coordCTD[:,0]+coordCTD[:,1]/60.+coordCTD[:,2]/3600.
lonCTD=-1*(coordCTD[:,3]+coordCTD[:,4]/60.+coordCTD[:,5]/3600.)

latCTD=np.unique(latCTD)

# Load data from temperature
nc=netcdf.Dataset(resdir+resfiles[0])
latitude = nc.variables['x'][:]
depth = -1*nc.variables['y'][:]
field = nc.variables['analyzed_field'][:]
nc.close()

# Load data from salinity
nc=netcdf.Dataset(resdir+resfiles[1])
latitude = nc.variables['x'][:]
depth = -1*nc.variables['y'][:]
field2 = nc.variables['analyzed_field'][:]
nc.close()

field=np.ma.masked_where(field==valex,field)
field2=np.ma.masked_where(field2==valex,field2)

# --------------
# Make the plot
# --------------

fig=plt.figure(figsize=(8,14))

# 1. Temperature
# ---------------

ax = host_subplot(221)
contour=ax.contour(latitude,depth,field,levels2plot,colors='k')
ax.clabel(contour,levels2plot[::2],inline=1,fmt='%1.f',fontsize=14,color='k')
ax.set_ylim((-depthmin,20))
ax.set_yticks(depth2plot)
ax.set_yticklabels(depth2plot2)
ax.set_xticks(np.arange(30.2,31.3,0.2))
ax.set_xticklabels('')
ax.set_ylabel('Depth (m)')

ax2 = ax.twin() # now, ax2 is responsible for "top" axis and "right" axis
ax2.set_xticks(sorted(latCTD))
ax2.set_xticklabels(stationlabels)
ax2.axis["right"].major_ticklabels.set_visible(False)
ax2.set_xlabel('Station no.')

# 2. Salinity
# ---------------

ax = host_subplot(222)
contour=ax.contour(latitude,depth,field2,levels2plot2,colors='k')
ax.clabel(contour,levels2plot2[::2],inline=1,fmt='%1.1f',fontsize=14,color='k')
ax.set_ylim((-depthmin,20))
ax.set_yticks(depth2plot)
ax.set_yticklabels('')
ax.set_xticks(np.arange(30.2,31.3,0.2))
ax.set_xticklabels('')


ax2 = ax.twin() # now, ax2 is responsible for "top" axis and "right" axis
ax2.set_xticks(sorted(latCTD))
ax2.set_xticklabels(stationlabels)
ax2.axis["right"].major_ticklabels.set_visible(False)
ax2.set_xlabel('Station no.')

# 3. Velocity
# ---------------

ax = host_subplot(223)
contour=ax.contour(latitude,depth,field,levels2plot,colors='k')
ax.clabel(contour,levels2plot[::2],inline=1,fmt='%1.f',fontsize=14,color='k')
ax.set_ylim((-depthmin,20))
ax.set_yticks(depth2plot)
ax.set_yticklabels(depth2plot2)
ax.set_xticks(np.arange(30.2,31.3,0.2))
ax.set_xlabel('Latitude ($^{\circ}$N)')
ax.set_ylabel('Depth (m)')

ax2 = ax.twin() # now, ax2 is responsible for "top" axis and "right" axis
ax2.set_xticks(sorted(latCTD))
ax2.set_xticklabels('')
ax2.axis["right"].major_ticklabels.set_visible(False)


# 4. Salinity
# ---------------

ax = host_subplot(224)
contour=ax.contour(latitude,depth,field2,levels2plot2,colors='k')
ax.clabel(contour,levels2plot2[::2],inline=1,fmt='%1.1f',fontsize=14,color='k')
ax.set_ylim((-depthmin,20))
ax.set_yticks(depth2plot)
ax.set_yticklabels('')
ax.set_xticks(np.arange(30.2,31.3,0.2))
ax.set_xlabel('Latitude ($^{\circ}$N)')

ax2 = ax.twin() # now, ax2 is responsible for "top" axis and "right" axis
ax2.set_xticks(sorted(latCTD))
ax2.set_xticklabels('')
ax2.axis["right"].major_ticklabels.set_visible(False)










plt.show()
plt.close()
