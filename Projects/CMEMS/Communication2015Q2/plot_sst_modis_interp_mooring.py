__author__ = 'ctroupin'


import glob
import os
import netCDF4
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from scipy import interpolate

figdir = "/data_local/Satellite/MODIS/data/daily/sst/terra/4mu/4km/figures/"
datadir = "/data_local/Satellite/MODIS/data/daily/sst/terra/4mu/4km/"
mooringfile = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/history/mooring/2plot/IR_TS_MO_61198.nc"

cmap = plt.cm.hot_r
fmin, fmax = 22., 29.
norm = colors.Normalize(vmin=fmin, vmax=fmax)
bounds = np.arange(fmin, fmax, 2.0)
coordinates = np.array((-6.75, 2., 35, 39.))
dlon, dlat = 2., 1.

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
            urcrnrlat=coordinates[3], lat_ts=0.5*(coordinates[2]+coordinates[3]), resolution='h')

with netCDF4.Dataset(mooringfile) as nc:
    latmooring = nc['LATITUDE'][:].mean()
    lonmooring = nc['LONGITUDE'][:].mean()

lonmooring2, latmooring2 = m(lonmooring, latmooring)
filelist = sorted(glob.glob(datadir+'*nc'))

for datafiles in filelist:
    ff = os.path.basename(datafiles)
    ff = ff.split('.')[0]

    with netCDF4.Dataset(datafiles) as nc:
        sst = nc.variables['sst4'][:]
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        timesst = nc.time_coverage_start

    # Interpolate the sst field onto the mooring position
    llon, llat = np.meshgrid(lon, lat)

    f = interpolate.interp2d(lon, lat, sst, kind='linear')
    sst_interp = f(lonmooring, latmooring)


    lon, lat = m(llon, llat)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m.pcolormesh(lon, lat, sst, cmap=cmap, norm=norm)
    if sst_interp<35:
        plt.text(lonmooring2, latmooring2, '  ' + str(np.round(sst_interp[0], 2)), ha='left', va='center')
        scat = m.scatter(lonmooring2, latmooring2, s=45, c=sst_interp, cmap=cmap, norm=norm)

    m.drawcoastlines(ax=ax, linewidth=0.2)
    m.fillcontinents(color = 'gray')
    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)
    cbar = plt.colorbar(scat, extend='both', shrink=0.7)
    cbar.set_label('$^{\circ}$C', rotation=0, horizontalalignment='left')
    cbar.set_ticks(bounds)
    plt.title(timesst[:10])
#    plt.savefig(figdir+ff)
    plt.show()
    plt.close()