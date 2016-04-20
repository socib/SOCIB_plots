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
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
import pysocibclient
import os
from osgeo import gdal

doplot = 1  # switch for the plot
cmap = plt.cm.RdYlBu_r

res = 'l'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'  # Coastline
coastfile = 'coatline_cartex.dat'
topodir = '/home/ctroupin/DataOceano/Bathymetry/'
topofile = 'medsea_bathymetry_etopo2.nc'
figdir = '/home/ctroupin/IMEDEA/Cartex2014/figures/deployment_time8/'  # where plots are saved
altimetryfile = "http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/madt/L4/2014/06/nrt_med_allsat_madt_uv_20140625_20140701.nc.gz"
altimetryfile2 = ("http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/" +
                 "madt/L4/2014/06/nrt_med_allsat_madt_h_20140625_20140701.nc.gz")

rcParams['contour.negative_linestyle'] = 'solid'

projectname = 'PERSEUS'

levels2plot = np.arange(-0.3, 0.30001, 0.02)


# region of interest
coordinates = np.array((-2, 6., 35.5, 41.001))
dlon, dlat = 1.0, 1.0

valex = -999.
tt = datetime.datetime(2014, 5, 25, 15, 0, 0)
tt_end = datetime.datetime(2014, 7, 31, 15, 0, 0)

time_min = int(tt.strftime('%s'))
time_max = int(tt_end.strftime('%s'))
delta_time = 10800.
time4interp = np.arange(time_min, time_max, delta_time)
ntimes = len(time4interp)

perseus_init = '2014-05-25T000000'

visiblefile = '/data_local/Satellite/Visible/nasa-worldview-2014-06-14.tiff'


gtif = gdal.Open(visiblefile)
gtif.GetProjectionRef()

# Set the plot axis limits to the proper map coordinates.
arr = gtif.ReadAsArray()
trans = gtif.GetGeoTransform()
extent = (trans[0], trans[0] + gtif.RasterXSize*trans[1],
          trans[3] + gtif.RasterYSize*trans[5], trans[3])

# prepare colormaps
vmin, vmax = 19, 25.5
sstmin, sstmax = 19, 25.5
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 1.0)
# Colormap
cmapsst = plt.cm.RdYlBu_r


# Load altimetry
with netcdf.Dataset(altimetryfile) as nc:
    lon_alti = nc.variables['lon'][:] - 360.
    lat_alti = nc.variables['lat'][:]
    u = nc.variables['u'][:].squeeze()
    v = nc.variables['v'][:].squeeze()

with netcdf.Dataset(altimetryfile2) as nc:
    lon_alti2 = nc.variables['lon'][:] - 360.
    lat_alti2 = nc.variables['lat'][:]
    adt = nc.variables['adt'][:].squeeze()

normalti = colors.Normalize(vmin=-0.15, vmax=0.15)

# Generate lists of platforms
socib_api = pysocibclient.ApiClient()
drifterlist = socib_api.list_platforms(init_datetime=perseus_init, instrument_type="surface_drifter")

loncoast, latcoast = np.loadtxt(coastdir + coastfile, usecols=(0, 1), unpack=True)
loncoast[loncoast == valex] = np.nan
latcoast[latcoast == valex] = np.nan

# Prepare plot
# ------------

if doplot:

    # Create directory
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    xc, yc = loncoast, latcoast

    # Where the text will be added
    xtext, ytext = -1.15, 39.5

ndrifter = len(drifterlist)
print 'Working on ' + str(ndrifter) + ' drifters'
print ' '

# Allocate
londrifter_interp = valex * np.ma.ones((ndrifter, ntimes))
latdrifter_interp = valex * np.ma.ones((ndrifter, ntimes))
WTEM_interp = valex * np.ma.ones((ndrifter, ntimes))
kk = 0

for drifter in drifterlist:

    drifter_opendap = drifter.product_list[-1].opendap
    nc = netcdf.Dataset(drifter_opendap)
    if nc.project == projectname:
        # print 'kk=',
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        gtime = nc.variables['time'][:]
        try:
            WTEM = nc.variables['WTEM'][:]
            WTEM_interp[kk, :] = np.interp(time4interp, gtime, WTEM)
        except KeyError:
            # print 'No variable WTEM'
            WTEM_interp[kk, :] = valex
        nc.close()

        londrifter_interp[kk, :] = np.interp(time4interp, gtime, lon)
        latdrifter_interp[kk, :] = np.interp(time4interp, gtime, lat)

        kk += 1

londrifter_interp = londrifter_interp[:kk, :]
latdrifter_interp = latdrifter_interp[:kk, :]

# Mask bad values
londrifter_interp = np.ma.masked_outside(londrifter_interp, coordinates[0]+1.0, coordinates[1])
latdrifter_interp = np.ma.masked_outside(latdrifter_interp, coordinates[2], coordinates[3])
latdrifter_interp = np.ma.masked_where(np.isnan(latdrifter_interp), latdrifter_interp)
londrifter_interp = np.ma.masked_where(np.isnan(londrifter_interp), londrifter_interp)

WTEM_interp = np.ma.masked_less(WTEM_interp, 10.)
WTEM_interp = np.ma.masked_where(np.isnan(WTEM_interp), WTEM_interp)



fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(coordinates[0], coordinates[1])
ax.set_ylim(coordinates[2], coordinates[3])
# m.drawparallels(np.arange(36., 41.0001, dlat), linewidth=0.,
#                 labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
# m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
#                 labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

plt.plot(xc, yc, 'k', lw=0.2, zorder=5)

p1 = plt.scatter(0.0, 0.0, c=15, s=0.0, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst)

speed = np.sqrt(u.data * u.data + v.data * v.data)
speed = np.ma.masked_greater(speed, 1.5)
plt.streamplot(lon_alti, lat_alti, u, v, color=adt, norm=normalti, arrowsize=2,
               density=12, linewidth=.25, cmap=plt.cm.RdBu_r, zorder=5)

plt.imshow(arr[:3,:,:].transpose((1, 2, 0)), extent=extent, zorder=3,alpha=1)

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

cbar = fig.colorbar(p1, cmap=cmapsst, norm=normsst, orientation='vertical', pad=0.025, aspect=15, shrink=0.75,
                    extend='both')
cbar.set_ticks(boundsst)
cbar.set_label(r'$^{\circ}$C', fontsize=16, rotation=0, ha='left')



f = 0
for i in range(24, ntimes):

    ii = str(i).zfill(3)
    print ii

    # plt.contour(lon_alti2, lat_alti2, adt, levels2plot, colors="w", linewidths=0.2, zorder=5, alpha=0.5)

    for n1 in range(0, len(londrifter_interp)):
        plt.plot(londrifter_interp[n1, i], latdrifter_interp[n1, i], 'ko', ms=0.1)

        try:
            plt.scatter(londrifter_interp[n1, 0:i-1], latdrifter_interp[n1, 0:i-1], c=WTEM_interp[n1, 0:i-1],
                      s=2, edgecolor='none', cmap=cmapsst, norm=normsst, zorder=6, alpha=0.75)
            plt.scatter(londrifter_interp[n1, i], latdrifter_interp[n1, i], c=WTEM_interp[n1, i],
                      s=5, edgecolor='none', cmap=cmapsst, norm=normsst, zorder=7, alpha=1)
        except AttributeError:
            f+=1
            # print 'No temp '

    texttime = time.strftime("%Y\n%b %d\n%H:%M", time.gmtime(time4interp[i]))
    props = dict(boxstyle='square', facecolor='white', alpha=0.25)
    thetext = ax.text(xtext, ytext, texttime, fontsize=20, zorder=7, ha='center', bbox=props)

    plt.savefig(figdir+ii, dpi=300, transparent=False, bbox_inches='tight')
    thetext.remove()

    # plt.show()

