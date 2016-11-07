#!/usr/bin/env python

import os
import glob
import numpy as np
import logging
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
import datetime, time, calendar
import matplotlib.text as text


from alborex_functions import *
import pysocibclient

logging.basicConfig(level=logging.INFO)

dlon, dlat =  1, 1.

# coordinates = [-1.5, 6., 36., 39.]
coordinates = [-1.25, 1.25, 36., 38.25]   # zoomed region

res = 'l'

# File and directory names
coastdir = '/home/ctroupin/PycharmProjects/SOCIB_plots/Missions/Alborex/'
coastfile = 'coastline_cartex.dat'
figdir = '/home/ctroupin/Publis/201502_Alborex/figures/Drifters/'
figname = 'Alborex_drifters_sst_zoom'


valex=-999

# Time for the interpolation and the drifters
#time_init='2014-01-01T000000'
time_init='2014-05-25T000000'
time_end='2014-07-15T000000'
tt = datetime.datetime( 2014, 5, 25, 1, 0, 0 )
tt2 = datetime.datetime( 2014, 5, 26, 1, 0, 0 )
tt_end = datetime.datetime( 2014, 6, 25, 15, 0, 0 )
time_min=int(tt.strftime('%s') )
time_min2=int(tt2.strftime('%s') )
time_max=int(tt_end.strftime('%s'))
delta_time=21600.
time4interp=np.arange(time_min, time_max, delta_time)
ntimes=len(time4interp)

# Colormap
cmapsst=plt.cm.RdYlBu_r

# Compute min and max values
sstmin,sstmax = 18., 24.
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1, 1.0)

logging.info("Loading coast")
loncoast,latcoast = alborex_load_coast(coastdir, coastfile, valex)

# Load drifters
# Generate lists of platforms
socib_api = pysocibclient.ApiClient()
drifterlist = socib_api.list_platforms(init_datetime=time_init,
                                       end_datetime=time_end,
                                       instrument_type="surface_drifter")

# Prepare the plot
fig,m, ax = prepare_map(coordinates,res)
loncoast,latcoast = m(loncoast,latcoast)

m.plot(loncoast, latcoast, 'k-', lw=0.25)

k=0
logging.info('Working on {0} drifters'.format(len(drifterlist)))
logging.info('Plotting drifter trajectories')
for drifter in drifterlist:

    logging.info('----------------')
    logging.info("Working on file {0}".format(drifter))
    drifter_opendap = drifter.product_list[-1].opendap
    logging.info(drifter_opendap)

    with netcdf.Dataset(drifter_opendap) as nc:

        # Check if drifter corresponds to Alborex
        if nc.project == 'PERSEUS':
            logging.info('Alborex data')

            lon, lat = nc.variables['LON'][:], nc.variables['LAT'][:]
            gtime = nc.variables['time'][:]
            # Select only good time
            if gtime[0] <= time_min2:
                goodtime = np.where((gtime<time_max)&(gtime>=time_min))
                lon,lat = lon[goodtime],lat[goodtime]
                londrifter, latdrifter=m(lon,lat)
                try:
                    WTEM = nc.variables['WTEM'][goodtime]
                    m.scatter(londrifter, latdrifter, c=WTEM, s=1.,
                              edgecolor='none', zorder=5, cmap=cmapsst, norm=normsst)
                    if (len(WTEM)>0):
                        pp=m.scatter(londrifter[1], latdrifter[1], c=WTEM[1],
                                     s=1, edgecolor='none', zorder=5,
                                     cmap=cmapsst, norm=normsst)
                except KeyError:
                    logging.info('No variable WTEM')
                    m.plot(londrifter, latdrifter,'ko', ms=0.25, zorder=4)
                if (len(londrifter)>0):
                    k=k+1
                    #m.plot(londrifter[0], latdrifter[0], 'ko', ms=3, zorder=6)
                    # m.plot(londrifter[-1], latdrifter[-1], 'ks', ms=3, zorder=6)
        else:
            logging.info('Not Alborex data')




cbar=fig.colorbar(pp,cmap=cmapsst, norm=normsst, orientation='vertical',
                  pad=0.025, aspect=15, shrink=0.6, extend='both')

cbar.set_ticks(boundsst)
cbar.set_label(r'$^{\circ}$C', fontsize=14, rotation=0, ha='left')

m.drawparallels(np.arange(coordinates[2],coordinates[3],dlat), linewidth=0.5,
                        labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=1)
m.drawmeridians(np.arange(coordinates[0],coordinates[1],dlon), linewidth=0.5,
                        labels=[0, 0, 0, 1], fontsize=16, zorder=1)

m.fillcontinents(ax=ax,color='w', zorder=2)
plt.title('2014-05-25 $-$ 2014-06-25', fontsize=20)


plt.savefig(os.path.join(figdir, figname), dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)

# plt.show()
plt.close()
