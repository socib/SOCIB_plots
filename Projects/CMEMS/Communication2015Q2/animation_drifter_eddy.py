#!/usr/bin/env python
#
# plot_deployments_time.py
#
# Plot the trajectories of the deployments:
# Interpolate positions every hour for each deployment
# Then plot
#
# ctroupin, June 2014
# --------------------------------------------------------------------------------

import numpy as np
import netCDF4 as netcdf
import time, datetime
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
import os

doplot = 1  # switch for the plot
cmap = plt.cm.RdYlBu_r

res = 'h'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'  # Coastline
coastfile = 'coatline_cartex.dat'
figdir = '/home/ctroupin/IMEDEA/Cartex2014/figures/EddyAnimation2/'  # where plots are saved

# region of interest
coordinates = np.array((3.5, 8.001, 36.5, 39))
dlon, dlat = 1.0, 1.0

valex = -999.
tt = datetime.datetime(2015, 6, 17, 11, 0, 0)
tt_end = datetime.datetime(2015, 8, 1, 2, 0, 0)

time_min = int(tt.strftime('%s'))
time_max = int(tt_end.strftime('%s'))
delta_time = 7200.
time4interp = np.arange(time_min, time_max, delta_time)
ntimes = len(time4interp)

tt1 = datetime.datetime(2015, 7, 1, 0, 0, 0)
tt2 = datetime.datetime(2015, 8, 1, 0, 0, 0)
tt1 = int(tt1.strftime('%s'))
tt2 = int(tt2.strftime('%s'))

# prepare colormaps
sstmin, sstmax = 23., 27.5
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 1.0)
# Colormap
cmapsst = plt.cm.RdYlBu_r


# Generate lists of platforms
drifterfile = "http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp055-ime_svp020/L1/2015/dep0002_drifter-svp055_ime-svp020_L1_2015-04-25.nc"

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
xtext, ytext = m(4.5, 38.5)


with netcdf.Dataset(drifterfile) as nc:
    lon = nc.variables['LON'][:]
    lat = nc.variables['LAT'][:]
    gtime = nc.variables['time'][:]
    try:
        WTEM = nc.variables['WTEM'][:]
        WTEM_interp = np.interp(time4interp, gtime, WTEM)
    except KeyError:
        print 'No variable WTEM'
        WTEM_interp = valex

londrifter_interp = np.interp(time4interp, gtime, lon)
latdrifter_interp = np.interp(time4interp, gtime, lat)

# Mask bad values
londrifter_interp = np.ma.masked_where(np.logical_or(time4interp >= time_max, time4interp <= time_min),
                                       londrifter_interp)
latdrifter_interp = np.ma.masked_where(np.logical_or(time4interp >= time_max, time4interp <= time_min),
                                       latdrifter_interp)
londrifter_interp = np.ma.masked_outside(londrifter_interp, coordinates[0], coordinates[1])
latdrifter_interp = np.ma.masked_outside(latdrifter_interp, coordinates[2], coordinates[3])


WTEM_interp = np.ma.masked_less(WTEM_interp, 10.)
WTEM_interp = np.ma.masked_where(np.isnan(WTEM_interp), WTEM_interp)

londrifter_interp, latdrifter_interp = m(londrifter_interp, latdrifter_interp)

def prepare_plot():
    fig = plt.figure(figsize=(8, 4.5))
    ax = fig.add_subplot(111)
    ax.set_xlim(coordinates[0], coordinates[1])
    ax.set_ylim(coordinates[2], coordinates[3])
    m.drawparallels(np.arange(36., 40., dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)


    m.plot(xc, yc, 'k', lw=0.2, zorder=3)

    p1 = m.scatter(0.0, 0.0, c=15,
                   s=0.0, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst)

    cbar = fig.colorbar(p1, cmap=cmapsst, norm=normsst, orientation='vertical', pad=0.025, aspect=15, shrink=0.99,
                        extend='both')
    cbar.set_ticks(boundsst)
    cbar.set_label(r'$^{\circ}$C', fontsize=16, rotation=0, ha='left')

    m.fillcontinents(ax=ax, color='grey', zorder=2)

    return fig


tempmin, tempmax = WTEM_interp.min(), WTEM_interp.max()
# plt.show()
f = 0
for i in range(359, ntimes):

    ii = str(i).zfill(3)
    print ii

    fig = prepare_plot()

    m.plot(londrifter_interp[i], latdrifter_interp[i], 'ko', ms=0.1)

    try:
        m.scatter(londrifter_interp[0:i], latdrifter_interp[0:i], c=WTEM_interp[0:i],
                  s=4, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst, alpha=1)
        m.scatter(londrifter_interp[i], latdrifter_interp[i], c=WTEM_interp[i],
                  s=7.5, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst, alpha=1)
    except AttributeError:
        f+=1
        # print 'No temp '

    # texttime = time.strftime("%Y\n%b %d\n%H:%M", time.gmtime(time4interp[i]))
    texttime = time.strftime("%Y-%m-%d %H:%M", time.gmtime(time4interp[i]))
    plt.title(texttime)

    a = plt.axes([0.175, 0.6, .25, .25])

    plt.scatter(time4interp[:i], WTEM_interp[:i], s=3, c= WTEM_interp[:i],
                cmap=cmapsst, norm=normsst, edgecolor='None')
    plt.setp(a, xlim=(time4interp[0],time4interp[-1]), ylim=(tempmin, tempmax), xticks=[tt1, tt2],
             xticklabels=['June', 'July'], yticks=[23., 25., 27.])



    # plt.text(xtext, ytext, texttime, fontsize=20, zorder=6, ha='center')

    plt.savefig(figdir+ii, dpi=300, transparent=False, bbox_inches='tight')
    # plt.show()
    plt.close(fig)
    plt.cla()
