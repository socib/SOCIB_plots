
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import glob
import os
from matplotlib  import colors
from alborex_functions import prepare_map

datadir = "/data_local/Satellite/MODIS/data/L2/CentralAmerica/"
figdir = "/data_local/Satellite/MODIS/figures/Papagayo/"
filelist = sorted(glob.glob(datadir + '*OC.nc'))
coordinates = [-103, -97, 10, 15.5]
dlon, dlat = 1., 1.
res = 'i'
# Colormap
cmapchloro=plt.cm.YlGnBu_r

# Compute min and max values
chloromin,chloromax = 0.05,0.6
normchloro= colors.Normalize(vmin=chloromin,vmax=chloromax)
boundchloro = np.arange(chloromin,chloromax+.001,1)
newticks = np.arange(0.0,0.6001,0.1)



for datafiles in filelist:
    print 'Workin on %s', datafiles


    # Load file content
    with netCDF4.Dataset(datafiles, 'r') as nc:
        lon = nc.groups['navigation_data'].variables['longitude'][:]
        lat = nc.groups['navigation_data'].variables['latitude'][:]
        chloro = nc.groups['geophysical_data'].variables['chlor_a'][:]
        qual = nc.groups['geophysical_data'].variables['l2_flags'][:]
        timechla = nc.time_coverage_start[:10]

    chloro = np.ma.masked_where(qual>1, chloro)

    # Plot
    fig, m, ax = prepare_map(coordinates,res)

    x,y = m(lon,lat)
    sstpcm = m.pcolormesh(x, y, chloro, cmap=cmapchloro, norm=normchloro)

    cbar = fig.colorbar(sstpcm, cmap=cmapchloro, norm=normchloro, orientation='vertical', pad=0.025,
                      aspect=15,shrink=1,extend='both')
    cbar.set_ticks(newticks)

    m.drawparallels(np.arange(coordinates[2] ,coordinates[3], dlat), linewidth=0.,
                            labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=1)
    m.drawmeridians(np.arange(coordinates[0], coordinates[1], dlon), linewidth=0.,
                            labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16,zorder=1)

    m.fillcontinents(ax=ax, color='0.5',zorder=2)
    m.drawcoastlines(ax=ax, linewidth=0.2)

    plt.title(timechla,fontsize=20)


    plt.savefig(figdir + 'chloro_' + timechla, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)


    # plt.show()
    plt.close()