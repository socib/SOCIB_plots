#!/usr/bin/python
__author__ = 'ctroupin'

import glob
import cf
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import shutil
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from matplotlib  import colors

doplot = 1

year = 2015
month = 7


figdir = "/home/ctroupin/DataOceano/MyOcean/figures/drifters/"
basedir = "/home/ctroupin/DataOceano/MyOcean/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/" + str(year) + str(month).zfill(2) + '/'
landfile = "/data_local/Bathymetry/ne_10m_land"

sstmin, sstmax = 5., 30.
cmapsst = plt.cm.RdYlBu_r
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1,2.5)
level2plotsst = np.arange(10, sstmax+.1,0.01)

def haversine(point1, point2):
    """Gives the distance between two points on earth.
    """
    earth_radius_km =  6371.
    lat1, lon1 = (np.deg2rad(coord) for coord in point1)
    lat2, lon2 = (np.deg2rad(coord) for coord in point2)
    dlat, dlon = (lat2 - lat1, lon2 - lon1)
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    great_circle_distance = 2 * np.arcsin(min(1,np.sqrt(a)))
    d = earth_radius_km * great_circle_distance
    return d

yearlist = (2014, 2015)
monthlist = range(1, 13)

m = Basemap(projection='robin',lon_0=0,resolution='c')

if doplot:
    fig = plt.figure(figsize=(12, 8))
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.00)

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
k, j = 0, 0
filelist = sorted(glob.glob(basedir+'*.nc'))
for datafiles in filelist:

    # print datafiles
    with netCDF4.Dataset(datafiles) as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]
        depth = nc.variables['DEPH'][:]
        distance = haversine( (lat[0], lon[0]), (lat[-1], lon[-1]))

        #print depth.shape[0]
        #print 'Distance = ' + str(distance)
        lon, lat = m(lon, lat)


        try:
            temperature = nc.variables['TEMP'][:,0]

            # print j
            if (doplot == 1) & (depth.shape[1] == 1):

                scat = m.scatter(lon.mean(), lat.mean(), s=distance/15., c=temperature.mean(), edgecolor='None', cmap=cmapsst, norm=normsst)


        except KeyError:
            # print 'No variable temperature in this file'
            #temperature = np.nan*np.ones_like(lat)
            k+=1
    j+=1

if doplot:

    cbar = plt.colorbar(scat, cmap=cmapsst, norm=normsst, extend='both', shrink=0.7)
    cbar.set_label('$^{\circ}$C', rotation=0, ha='left')
    plt.title('Temperature from surface drifters\n' + str(year) + '-' + str(month).zfill(2))
    #plt.savefig(figdir + 'drifter_scatter')
    plt.show()
    #plt.close()

print 'Number of files: ' + str(len(filelist))
print 'Number of files without temperature: ' + str(k)





