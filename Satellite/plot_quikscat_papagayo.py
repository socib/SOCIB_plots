#!/usr/bin/env python
#
# download_plot_windsat.py
#
# Download is done with script download_wind_podaac
# ---------------------------------------------------------------------------------------

import os
import glob
import shutil
import numpy as np
from mpl_toolkits.basemap import Basemap
from wind_plot_functions import *
from scipy.interpolate import griddata

databasedir = '/data_local/Satellite/Wind/QuikScat/L2/'
datadir2 = databasedir + '2keep/'
figdir = '/data_local/Satellite/Wind/QuikScat/figures/papagayo/'

coordinates = [-110, -92, 6, 20.001, 2.0, 2.]

if not os.path.exists(datadir2):
    os.makedirs(datadir2)

if not os.path.exists(figdir):
    os.makedirs(figdir)

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2] + 0.1,
            urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
            lat_ts=0.5 * (coordinates[2] + coordinates[3]),
            resolution='l')

# Create grid for interpolation

lon2interp, lat2interp = np.meshgrid(np.arange(coordinates[0], coordinates[1], 1./8.),
                                     np.arange(coordinates[2], coordinates[3], 1./8.))

filelist = sorted(glob.glob(databasedir + 'qs_l2b_*nc'))
for datafiles in filelist:

    # Load file content
    lon, lat, uwind, vwind, windtime = read_L2_wind_quikscat(datafiles, coordinates)

    # Make the plot if data in study region
    if len(uwind)>0:
        shutil.copy2(datafiles, datadir2)

filelist2 = sorted(glob.glob(datadir2 + 'qs_l2b_*nc'))
for datafiles in filelist2[0:2]:
    print datafiles

    lon, lat, uwind, vwind, windtime = read_L2_wind_quikscat(datafiles, coordinates)

    # interpolate on regular grid! (use griddata)
    uwind_interp = griddata((lon,lat),uwind,(lon2interp,lat2interp))
    vwind_interp = griddata((lon,lat),vwind,(lon2interp,lat2interp))

    plot_wind_sat_interp(datafiles, coordinates, lon2interp, lat2interp, uwind_interp, vwind_interp, windtime, figdir, m)



