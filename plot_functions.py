import os
import numpy as np
import netCDF4 as netcdf
import cf
import matplotlib.pyplot as plt
import seawater
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_profiler_TS(temp, psal, depth, figname=None, **kwargs):
    fig = plt.figure(figsize=(10, 10))
    scat = plt.scatter(psal, temp, s=5, c=depth, edgecolor='None', **kwargs)
    cb = plt.colorbar(extend='max')
    # Compute density for T-S couples
    density2plot = np.arange(24, 30, .5)
    tvec = np.arange(np.nanmin(temp), np.nanmax(temp), 0.1)
    svec = np.arange(np.nanmin(psal), np.nanmax(psal), 0.05)
    ssvec, ttvec = np.meshgrid(svec, tvec)
    density = seawater.eos80.dens0(ssvec,ttvec) - 1000.0
    cont = plt.contour(svec, tvec, density, levels=density2plot, colors='.65',
                       linestyles='dashed', lineswidth=1.0)
    plt.clabel(cont, inline=True, fmt='%1.1f')
    
    plt.xlabel('Salinity')
    plt.ylabel('Temperature\n($^{\circ}$C)', rotation=0, ha='right')
    vmin, vmax = scat.get_clim()
    cb.set_ticks(np.linspace(vmin, vmax, 6))
    cb.set_label('Depth\n(m)', rotation=0, ha='left')
    if figname is None:
        plt.show()
    else:
        plt.savefig(figname)
    plt.close()

def plot_field_map(m, lon, lat, field, coordinates, Nticks, figname=None, **kwargs):

    lon, lat = m(lon, lat)

    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pcm = m.pcolormesh(lon, lat, field, **kwargs)
    
    m.drawcoastlines(linewidth=.5, zorder=3)
    m.fillcontinents(zorder=2)
    m.drawmeridians(np.linspace(coordinates[0], coordinates[1], Nticks), labels=[False,False,False,True], zorder=1)
    m.drawparallels(np.linspace(coordinates[2], coordinates[3], Nticks), labels=[True,False,False,False], zorder=1)

    cbar = plt.colorbar(pcm, extend='both', cax=cax)
    cbar.set_label('^{\circ}C', rotation=0, horizontalalignement='left')

    if figname is None:
        plt.show()
    else:
        plt.savefig(figname)
    plt.close()
