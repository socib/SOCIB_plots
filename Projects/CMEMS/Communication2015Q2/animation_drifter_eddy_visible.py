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
from osgeo import gdal


doplot = 1  # switch for the plot
cmap = plt.cm.RdYlBu_r

res = 'h'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'  # Coastline
coastfile = 'coatline_cartex.dat'
figdir = '/home/ctroupin/IMEDEA/Cartex2014/figures/EddyAnimationVisible2/'  # where plots are saved
visiblefile = '/data_local/Satellite/Visible/nasa-worldview-2014-06-14b.tiff'


gtif = gdal.Open(visiblefile)
gtif.GetProjectionRef()

#Set the plot axis limits to the proper map coordinates.
arr = gtif.ReadAsArray()
trans = gtif.GetGeoTransform()
extent = (trans[0], trans[0] + gtif.RasterXSize*trans[1],
          trans[3] + gtif.RasterYSize*trans[5], trans[3])

print extent

# region of interest
coordinates = np.array((2, 7.5, 36.5, 40.3))
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
sstmin, sstmax = 23., 29.001
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 1.0)
# Colormap
cmapsst = plt.cm.RdYlBu_r


# Generate lists of platforms
drifterfile = "http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp055-ime_svp020/L1/2015/dep0002_drifter-svp055_ime-svp020_L1_2015-04-25.nc"
altimetryfile = "http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/madt/L4/2015/07/nrt_med_allsat_madt_uv_20150708_20150714.nc.gz"

loncoast, latcoast = np.loadtxt(coastdir + coastfile, usecols=(0, 1), unpack=True)
loncoast[loncoast == valex] = np.nan
latcoast[latcoast == valex] = np.nan

if not os.path.exists(figdir):
    os.makedirs(figdir)

xc, yc = loncoast, latcoast

with netcdf.Dataset(altimetryfile) as nc:
    lon_alti = nc.variables['lon'][:] - 360.
    lat_alti = nc.variables['lat'][:]
    u = nc.variables['u'][:].squeeze()
    v = nc.variables['v'][:].squeeze()


# Where the text will be added
xtext, ytext = 4.5, 38.5

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

tempmin, tempmax = WTEM_interp.min(), WTEM_interp.max()
# plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(coordinates[0], coordinates[1])
ax.set_ylim(coordinates[2], coordinates[3])

plt.tick_params(
axis='x',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='off',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelbottom='off') # labels along the bottom edge are off

plt.tick_params(
axis='y',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
left='off',      # ticks along the bottom edge are off
right='off',         # ticks along the top edge are off
labelleft='off') # labels along the bottom edge are off
# m.drawparallels(np.arange(36., 40., dlat), linewidth=0.,
#                 labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
# m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
#                 labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)
#

# plt.plot(xc, yc, 'k', lw=0.2, zorder=3)
plt.streamplot(lon_alti, lat_alti, u, v, color='0.85',
               arrowsize=2,
               density=10, linewidth=.25, cmap=plt.cm.gray_r, zorder=5)

plt.imshow(arr[:3,:,:].transpose((1, 2, 0)), extent=extent, zorder=3, alpha=1)
p1 = plt.scatter(0.0, 0.0, c=15,
                   s=0.0, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst)

cbar = fig.colorbar(p1, cmap=cmapsst, norm=normsst, orientation='vertical', pad=0.025, aspect=15, shrink=0.75,
                    extend='both')
cbar.set_ticks(boundsst)
cbar.set_label(r'$^{\circ}$C', fontsize=16, rotation=0, ha='left')


f = 0
for i in range(1, ntimes):

    ii = str(i).zfill(3)
    print ii

    ax = fig.add_subplot(111)

    plt.scatter(londrifter_interp[0:i], latdrifter_interp[0:i], c=WTEM_interp[0:i],
                s=4, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst, alpha=1)
    plt.scatter(londrifter_interp[i], latdrifter_interp[i], c=WTEM_interp[i],
                s=7.5, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst, alpha=1)
    f += 1

    a = plt.axes([0.175, 0.25, .25, .25])
    plt.scatter(time4interp[:i], WTEM_interp[:i], s=3, c= WTEM_interp[:i],
                cmap=cmapsst, norm=normsst, edgecolor='None')
    plt.setp(a, xlim=(time4interp[0],time4interp[-1]), ylim=(tempmin, tempmax), xticks=[tt1, tt2],
             xticklabels=['June', 'July'], yticks=[23., 25., 27.])
    a.tick_params(axis='x', colors='white')
    a.tick_params(axis='y', colors='white')

    # texttime = time.strftime("%Y\n%b %d\n%H:%M", time.gmtime(time4interp[i]))
    texttime = time.strftime("%Y-%m-%d %H:%M", time.gmtime(time4interp[i]))
    plt.title(texttime, color='w')
    # plt.text(xtext, ytext, texttime, fontsize=20, zorder=6, ha='center')

    plt.savefig(figdir+ii, dpi=300, transparent=False, bbox_inches='tight')
#     plt.show()
#    plt.close(fig)
#    plt.cla()
