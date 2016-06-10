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
import pysocibclient
import os

doplot = 1  # switch for the plot
cmap = plt.cm.RdYlBu_r

res = 'h'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'  # Coastline
coastfile = 'coatline_cartex.dat'
topodir = '/home/ctroupin/DataOceano/Bathymetry/'
topofile = 'medsea_bathymetry_etopo2.nc'
figdir = '/home/ctroupin/IMEDEA/Cartex2014/figures/deployment_time6/'  # where plots are saved

projectname = 'PERSEUS'

# region of interest
coordinates = np.array((-2, 6., 35.5, 41.001))
dlon, dlat = 1.0, 1.0

valex = -999.
tt = datetime.datetime(2014, 5, 25, 15, 0, 0)
tt_end = datetime.datetime(2014, 7, 31, 15, 0, 0)

time_min = int(tt.strftime('%s'))
time_max = int(tt_end.strftime('%s'))
delta_time = 7200.
time4interp = np.arange(time_min, time_max, delta_time)
ntimes = len(time4interp)

perseus_init = '2014-05-25T000000'


# prepare colormaps
vmin, vmax = 18., 24.5
sstmin, sstmax = 18., 24.5
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 1.0)
# Colormap
cmapsst = plt.cm.RdYlBu_r


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

    # Create basemap then used for the plot
    m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
                urcrnrlat=coordinates[3], lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)

    # lon2plot,lat2plot=np.arange(np.ceil(lonmin),lonmax+0.0001,dlon),np.arange(latmin,latmax+0.001,dlat),

    cmap = plt.cm.Spectral_r
    # cmap=plt.cm.gist_ncar

    xc, yc = m(loncoast, latcoast)

    # Where the text will be added
    xtext, ytext = m(-1.15, 40.)

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
        print 'kk=', kk
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        gtime = nc.variables['time'][:]
        try:
            WTEM = nc.variables['WTEM'][:]
            WTEM_interp[kk, :] = np.interp(time4interp, gtime, WTEM)
        except KeyError:
            print 'No variable WTEM'
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

londrifter_interp, latdrifter_interp = m(londrifter_interp, latdrifter_interp)

WTEM_interp = np.ma.masked_less(WTEM_interp, 10.)
WTEM_interp = np.ma.masked_where(np.isnan(WTEM_interp), WTEM_interp)


def prepare_plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(coordinates[0], coordinates[1])
    ax.set_ylim(coordinates[2], coordinates[3])
    m.drawparallels(np.arange(36., 41.0001, dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)


    m.plot(xc, yc, 'k', lw=0.2, zorder=5)

    p1 = m.scatter(0.0, 0.0, c=15,
                   s=0.0, edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst)

    cbar = fig.colorbar(p1, cmap=cmapsst, norm=normsst, orientation='vertical', pad=0.025, aspect=15, shrink=0.95,
                        extend='both')
    cbar.set_ticks(boundsst)
    cbar.set_label(r'$^{\circ}$C', fontsize=16, rotation=0, ha='left')

    m.fillcontinents(ax=ax, color='grey', zorder=4)

    return fig, ax, m, cbar



# plt.show()
f = 0
for i in range(685, ntimes):

    ii = str(i).zfill(3)
    print ii
    fig, ax, m, cbar = prepare_plot()


    for n1 in range(0, len(londrifter_interp)):
        m.plot(londrifter_interp[n1, i], latdrifter_interp[n1, i], 'ko', ms=0.1)

        try:
            m.scatter(londrifter_interp[n1, 0:i-1], latdrifter_interp[n1, 0:i-1], c=WTEM_interp[n1, 0:i-1],
                      s=4, edgecolor='none', cmap=cmapsst, norm=normsst, zorder=2)
            m.scatter(londrifter_interp[n1, i], latdrifter_interp[n1, i], c=WTEM_interp[n1, i],
                      s=7.5, edgecolor='none', cmap=cmapsst, norm=normsst, zorder=3)
        except AttributeError:
            f+=1
            # print 'No temp '


    texttime = time.strftime("%Y\n%b %d\n%H:%M", time.gmtime(time4interp[i]))
    plt.text(xtext, ytext, texttime, fontsize=20, zorder=6, ha='center')

    plt.savefig(figdir+ii, dpi=300, transparent=False, bbox_inches='tight')
    # plt.show()
    plt.close()
