# Plot the coastline corresponding to
# the selected region
#

import glob
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from osgeo import gdal
import netCDF4 as netcdf
from matplotlib  import colors
import matplotlib.patheffects as path_effects
gdal.UseExceptions()

# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=-6,42,30.,48.,0.2,0.2

# create figure
fig = plt.figure()
ax = plt.subplot(111)

# set up orthographic map projection with
map = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin,\
            urcrnrlon=lonmax,urcrnrlat=latmax,  \
            lat_ts=0.5*(latmin+latmax),\
            resolution='i')

map.fillcontinents()
plt.show()
plt.close()
