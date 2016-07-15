'''
plot_colorbars_alone.py

Make a colorbar (or legend) as a separate figure.
for temperature, salinity and velocity
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from matplotlib  import colors

figtype='.png'
resolution=200
figdir='/home/ctroupin/SOCIB/Facilities/Modelling/figures/colorbars/'
fs=24           # fontsize

# Min and max values for each variables
tmin,tmax,tmin2,tmax2,dt=10.,28.,283.15,301.15,3.
smin,smax,ds=36.5,38.5,0.5
vmin,vmax,dv=0.0,1.0,0.2
wmin,wmax,dw=0.5,3.,0.5
windmin,windmax,dwind=0.0,30.,5.0


# Color map
cmap=plt.cm.RdYlBu_r
cmap=plt.cm.hot_r

jet_vals = np.array([[ 0.13333,0.9509,0.13333,1. ],         # Forest Green
                     [ 1,1,0,1.],
                     [ 1,0.0,0.0,1.]])

cmap = colors.LinearSegmentedColormap.from_list("cmapredorangered", jet_vals)

# define function



# Temperature
fig = plt.figure(figsize=(8,3.5))
fig.patch.set_facecolor('white')
fig.patch.set_alpha(0.)

ax1 = fig.add_axes([0.05, 0.30, 0.9, 0.125])
norm = mpl.colors.Normalize(vmin=tmin,vmax=tmax)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal',extend='both')
cb1.set_label('$^{\circ}$C',fontsize=fs)
cb1.set_ticks(np.arange(tmin,tmax+0.0001,dt))
cb1.ax.tick_params(labelsize=fs)
#cb1.ax.xaxis.set_ticks_position('top')
cb1.ax.set_xlim(tmin, tmax)

ax2 = ax1.twiny()
cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm,orientation='horizontal',extend='both')
###cb2.ax.xaxis.set_ticks_position('bottom')
##cb2.set_label('K',fontsize=fs)
cb2.set_ticks(np.arange(tmin,tmax+0.0001,dt))
cb2.set_ticklabels(np.arange(tmin,tmax+0.0001,dt))
cb2.ax.tick_params(labelsize=fs)
cb2.ax.set_xlim(tmin2, tmax2)
##plt.savefig(figdir+'temp_colorbar3'+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)
plt.show()
plt.close()

### Salinity
##fig = plt.figure(figsize=(8,1.5))
##fig.patch.set_facecolor('white')
##fig.patch.set_alpha(0.)
##ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
##norm = mpl.colors.Normalize(vmin=smin,vmax=smax)
##cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
##                                   norm=norm,
##                                   orientation='horizontal',extend='both')
##cb1.set_ticks(np.arange(smin,smax+0.0001,ds))
##cb1.ax.tick_params(labelsize=fs) 
##plt.savefig(figdir+'salt_colorbar3'+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)
###plt.show()
##plt.close()
##
##
### Velocity
##fig = plt.figure(figsize=(8,1.5))
##fig.patch.set_facecolor('white')
##fig.patch.set_alpha(0.)
##ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
##norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
##cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
##                                   norm=norm,
##                                   orientation='horizontal',extend='max')
##cb1.set_label('m/s',fontsize=fs)
##cb1.set_ticks(np.arange(vmin,vmax+0.0001,dv))
##cb1.ax.tick_params(labelsize=fs) 
##plt.savefig(figdir+'vel_colorbar3'+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)
###plt.show()
##plt.close()
##
### Wave height
##fig = plt.figure(figsize=(8,1.5))
##fig.patch.set_facecolor('white')
##fig.patch.set_alpha(0.)
##ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
##norm = mpl.colors.Normalize(vmin=wmin,vmax=wmax)
##cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
##                                   norm=norm,
##                                   orientation='horizontal',extend='both')
##cb1.set_label('m',fontsize=fs)
##cb1.set_ticks(np.arange(wmin,wmax+0.0001,dw))
##cb1.ax.tick_params(labelsize=fs) 
##plt.savefig(figdir+'wave_sig_height_colorbar'+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)
###plt.show()
##plt.close()
##
##
### Wave height
##fig = plt.figure(figsize=(8,1.5))
##fig.patch.set_facecolor('white')
##fig.patch.set_alpha(0.)
##ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
##norm = mpl.colors.Normalize(vmin=windmin,vmax=windmax)
##cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
##                                   norm=norm,
##                                   orientation='horizontal',extend='max')
##cb1.set_label('km/h',fontsize=fs)
##cb1.set_ticks(np.arange(windmin,windmax+0.0001,dwind))
##cb1.ax.tick_params(labelsize=fs) 
##plt.savefig(figdir+'wind_colorbar3'+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)
###plt.show()
##plt.close()
##
### Wave period
##fig = plt.figure(figsize=(8,1.5))
##fig.patch.set_facecolor('white')
##fig.patch.set_alpha(0.)
##ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
##norm = mpl.colors.Normalize(vmin=windmin,vmax=windmax)
##cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
##                                   norm=norm,
##                                   orientation='horizontal',extend='max')
##cb1.set_label('$s$',fontsize=fs)
##cb1.set_ticks(np.arange(windmin,windmax+0.0001,dwind))
##cb1.ax.tick_params(labelsize=fs) 
##plt.savefig(figdir+'wind_colorbar3'+figtype, dpi=resolution, facecolor='w', edgecolor='w',
##      transparent=True, bbox_inches='tight', pad_inches=0.)
###plt.show()
##plt.close()
##
##
##
