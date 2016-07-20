#!/usr/bin/python

import os
import numpy as np
import netCDF4 as netcdf
import cf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_field_map(m, lon, lat, field, coordinates, Nticks, figname, cmap):

    lon, lat = m(lon, lat)

    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    divider = make_axes_locatable(ax)
    
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pcm = m.pcolormesh(lon, lat, field, vmin=23, vmax=27, cmap=cmap)
    
    m.drawcoastlines(linewidth=.5, zorder=3)
    m.fillcontinents(zorder=2)
    m.drawmeridians(np.linspace(coordinates[0], coordinates[1], Nticks), labels=[False,False,False,True], zorder=1)
    m.drawparallels(np.linspace(coordinates[2], coordinates[3], Nticks), labels=[True,False,False,False], zorder=1)

    cbar = plt.colorbar(pcm, extend='both', cax=cax)
    cbar.set_label('^{\circ}C', rotation=0, horizontalalignement='left')
    # plt.show()
    plt.savefig(figname)
    plt.close()
