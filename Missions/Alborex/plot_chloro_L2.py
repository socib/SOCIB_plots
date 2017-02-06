#!/usr/bin/env python
#
# -------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import datetime, time, calendar
import scipy.io
import matplotlib.text as text

from alborex_functions import *

# File and directory names

dlon,dlat =  1.0,1.0
coordinates = np.array((-3,1.,35.,38.2))
res = 'l'
coastdir='/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'
coastfile='coastline_cartex.dat'

bathydir = '/home/ctroupin/IMEDEA/Cartex2014/data/bathymetry/'
bathyfile = 'topo_gebco_medsea.nc'
ctdfile = ('http://thredds.socib.es/thredds/dodsC/research_vessel/ctd/'
           'socib_rv-scb_sbe9002/L1/2014/'
           'dep0007_socib-rv_scb-sbe9002_L1_2014-05-25.nc')
gliderfile1 = ('http://thredds.socib.es/thredds/dodsC/auv/glider/'
               'ideep00-ime_sldeep000/L2/2014/'
               'dep0012_ideep00_ime-sldeep000_L2_2014-05-25_data_dt.nc')
gliderfile2 = ('http://thredds.socib.es/thredds/dodsC/auv/glider/'
               'icoast00-ime_slcost000/L2/2014/'
               'dep0005_icoast00_ime-slcost000_L2_2014-05-25_data_dt.nc')

chlorodir = '/data_local/Satellite/MODIS/data/L2/Alborex/OC/'
figdir = '/home/ctroupin/Publis/201502_Alborex/figures/chloro/'

chlorofilelist = sorted(glob.glob(chlorodir + '*.nc'))
nfiles = len(chlorofilelist)
valex = 999

# Colormap
cmapchloro = plt.cm.YlGnBu_r

# Compute min and max values
chloromin,chloromax = 0.05,0.6
normchloro = colors.LogNorm(vmin=chloromin,vmax=chloromax)
boundchloro = np.arange(chloromin,chloromax+.001,1)

newticks = np.array((0.06, 0.1, 0.3, 1.0,))
newlabels = np.array((0.06, 0.1, 0.3, 1.0,))

### Load coast
loncoast,latcoast = alborex_load_coast(coastdir,coastfile,valex)

for f in range(0,1):

    fig, m, ax = prepare_map(coordinates, res)
    loncoast, latcoast = m(loncoast, latcoast)

    # Load L2 data

    with netcdf.Dataset(chlorofilelist[f], 'r') as nc:
      lon = nc.variables['Navigation_Data_longitude'][:]
      lat = nc.variables['Navigation_Data_latitude'][:]
      chla = nc.variables['Geophysical_Data_chlor_a'][:]
      timechla = nc.time_coverage_start

    # Build date
    datechla = timechla[:4] + '-' + timechla[4:6] + '-' + timechla[6:8]

    # Mask
    NN=1
    np.ma.masked_inside(chla,0,10)

    #x,y=m(lon[latstart:-latend:NN,lonstart:-lonend:NN],lat[latstart:-latend:NN,lonstart:-lonend:NN])
    x,y=m(lon,lat)
    chloropcm = m.pcolormesh(x, y, chla,
                             cmap=cmapchloro, norm=normchloro,
                             edgecolor='none')

    m.plot(lonCTD, latCTD,'ko',ms=3,zorder=2)
    m.plot(loncoast,latcoast,'r-',lw=0.5)

    cbar=fig.colorbar(chloropcm, cmap=cmapchloro ,norm=normchloro,
                      orientation='vertical', pad=0.025, aspect=15,
                      shrink=1, extend='max')

    cbar.set_ticks(newticks)
    cbar.set_ticklabels(newlabels)

    m.drawparallels(np.arange(coordinates[2],coordinates[3],dlat), linewidth=0.5,
                            labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=1)
    m.drawmeridians(np.arange(coordinates[0],coordinates[1],dlon), linewidth=0.5,
                            labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16,zorder=1)
    m.drawmapscale(.5,35.5,1,37,50, barstyle='simple', units='km', fontsize=12,zorder=3)

    m.fillcontinents(ax=ax,color='w',zorder=2)
    m.drawcoastlines(ax=ax)

    plt.title(datechla,fontsize=20)

    figname=os.path.basename(chlorofilelist[f])[:-3]
##    plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',
##                 transparent=False, bbox_inches='tight', pad_inches=0.1)


    plt.show()
    plt.close()
