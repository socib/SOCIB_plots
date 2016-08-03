#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from matplotlib import colors
from scipy.interpolate import griddata

doplot = 1

year = 2015
month = 7


figdir = "/home/ctroupin/DataOceano/MyOcean/figures/drifters/"
basedir = "/home/ctroupin/DataOceano/MyOcean/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/" + str(year) + str(month).zfill(2) + '/'
landfile = "/data_local/Bathymetry/ne_10m_land"

sstmin, sstmax = 5., 30.
cmapsst = plt.cm.RdYlBu_r
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1, 2.5)
level2plotsst = np.arange(10, sstmax+.1, 0.01)

lon_interp, lat_interp = np.meshgrid(np.arange(-180, 180.), np.arange(-80, 80.))


yearlist = (2014, 2015)
monthlist = range(1, 13)

if doplot:
    fig = plt.figure(figsize=(12, 8))
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.00)

    m = Basemap(projection='robin',lon_0=0,resolution='c')
    shp_info = m.readshapefile(landfile, 'scalerank', drawbounds=True)
    ax = plt.gca()
    ax.cla()

    paths = []
    for line in shp_info[4]._paths:
        paths.append(Path(line.vertices, codes=line.codes))

    coll = PathCollection(paths, linewidths=0, facecolors='grey', zorder=2)

    m = Basemap(projection='robin',lon_0=0,resolution='c')
    # drawing something seems necessary to 'initiate' the map properly
    m.drawcoastlines(color='white', zorder=0)
    ax = plt.gca()
    ax.add_collection(coll)

# Loop on files
k = 0
lon, lat, temperature = np.array([]), np.array([]), np.array([])

filelist = sorted(glob.glob(basedir+'*.nc'))
for datafiles in filelist:

    # print datafiles
    with netCDF4.Dataset(datafiles) as nc:
        lon = np.append(lon, nc.variables['LONGITUDE'][:])
        lat = np.append(lat, nc.variables['LATITUDE'][:])

#        lon2, lat2 = m(lon, lat)

        try:
            temperature = np.append(temperature, nc.variables['TEMP'][:,0])
        except KeyError:
            temperature = np.append(temperature, np.nan*np.ones_like(nc.variables['LATITUDE'][:]))
            # print 'No variable temperature in this file'
            # temperature = np.nan*np.ones_like(lat)
            k+=1

print lon.shape
print lat.shape
print temperature.shape

temperature_interp = griddata((lon, lat), temperature, (lon_interp, lat_interp), method='linear')
temperature_interp = np.ma.masked_where(np.isnan(temperature_interp), temperature_interp)
temperature = np.ma.masked_where(np.isnan(temperature), temperature)
temperature = np.ma.masked_outside(temperature, -2.0, 35.)
temperature_interp = np.ma.masked_outside(temperature_interp, -2.0, 35.)
print temperature_interp.min()
print temperature_interp.max()
print temperature.min()
print temperature.max()

# Loop on interpolation grid


lon_interp, lat_interp = m(lon_interp, lat_interp)
lon, lat = m(lon, lat)

m.pcolormesh(lon_interp, lat_interp, temperature_interp, cmap=cmapsst, norm=normsst)
scat = m.scatter(lon, lat, s=2, c=temperature, edgecolor='None', cmap=cmapsst, norm=normsst, zorder=4)

cbar = plt.colorbar(scat, cmap=cmapsst, norm=normsst, extend='both', shrink=0.7)
cbar.set_label('$^{\circ}$C', rotation=0, ha='left')
plt.title('Temperature from surface drifters\n' + str(year) + '-' + str(month).zfill(2))
plt.savefig(figdir + 'drifter_interp')
# plt.show()
plt.close()



print 'Number of files: ' + str(len(filelist))
print 'Number of files without temperature: ' + str(k)





