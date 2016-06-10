#!/usr/bin/env python
#
# ctroupin, August 2015
# --------------------------------------------------------------------------------

import numpy as np
import netCDF4 as netcdf
import time, datetime
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
import os
import scipy.io

doplot = 1  # switch for the plot
cmap = plt.cm.RdYlBu_r

res = 'h'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'  # Coastline
coastfile = 'coatline_cartex.dat'
figdir = '/home/ctroupin/IMEDEA/Cartex2014/figures/EddyAnimationOGS/'  # where plots are saved

# region of interest
coordinates = np.array((2, 6., 35.5, 41))
dlon, dlat = 1.0, 1.0

valex = -999.
tt = datetime.datetime(2014, 6, 24, 11, 0, 0)
tt_end = datetime.datetime(2015, 8, 3, 12, 0, 0)

time_min = int(tt.strftime('%s'))
time_max = int(tt_end.strftime('%s'))


# Generate lists of platforms
drifterfile = "/home/ctroupin/Projects/201501_InsTAC/DrifterOGS/136014.mat"



loncoast, latcoast = np.loadtxt(coastdir + coastfile, usecols=(0, 1), unpack=True)
loncoast[loncoast == valex] = np.nan
latcoast[latcoast == valex] = np.nan

if not os.path.exists(figdir):
    os.makedirs(figdir)

# Create basemap then used for the plot
m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
            urcrnrlat=coordinates[3], lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)

xc, yc = m(loncoast, latcoast)

# Where the text will be added
xtext, ytext = m(2.5, 38.5)

mat = scipy.io.loadmat(drifterfile)
lon = mat['Lon'][:]
lat = mat['Lat'][:]
gtime = mat['Time'][:]




londrifter_interp, latdrifter_interp = m(lon, lat)

def prepare_plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(coordinates[0], coordinates[1])
    ax.set_ylim(coordinates[2], coordinates[3])
    m.drawparallels(np.arange(36., 40., dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

    m.fillcontinents(ax=ax, color='grey', zorder=2)

    return fig



f = 0
for i in range(0, 1):

    ii = str(i).zfill(3)
    print ii

    fig = prepare_plot()

    m.plot(londrifter_interp[i], latdrifter_interp[i], 'ko', ms=10)

    texttime = time.strftime("%Y-%m-%d %H:%M", time.gmtime(gtime[i]))
    plt.title(texttime)




    # plt.text(xtext, ytext, texttime, fontsize=20, zorder=6, ha='center')

    plt.savefig(figdir+ii, dpi=300, transparent=False, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    plt.cla()
