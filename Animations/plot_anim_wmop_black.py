__author__ = 'ctroupin'

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import cf
import calendar
from matplotlib import colors
from matplotlib.path import Path
from matplotlib.collections import PathCollection

figdir = "/data_local/DataOceano/SOCIB/WMOP/figures/4animations/4/"


yearmin, yearmax = 2015, 2015
monthmin, monthmax = 1, 1
vmin, vmax, dvar = 36., 38.5, 0.5

variable = 'sea_surface_salinity'
cmap = plt.cm.RdYlBu_r
cmap.set_bad ('k',1.0)

coordinates = np.array((-5.5, 9.001, 34.75, 44.))
dlon, dlat = 3.0, 2.0

m = Basemap(llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
            urcrnrlon=coordinates[1], urcrnrlat=coordinates[3], resolution='l')

findex = 0

landfile = "/data_local/Bathymetry/ne_10m_land"
shp_info = m.readshapefile(landfile, 'scalerank', drawbounds=True)

def makeplot(field2plot, lon_ts, lat_ts, vmin, vmax, m, findex, shp_info):

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    bounds = np.arange(vmin, vmax + .0001, dvar)

    ax = plt.gca()
    ax.cla()

    paths = []
    for line in shp_info[4]._paths:
        paths.append(Path(line.vertices, codes=line.codes))
    coll = PathCollection(paths, linewidths=0, facecolors='black', zorder=2)
    ax.add_collection(coll)

    m.drawcoastlines(ax=ax, linewidth=0.05)
    m.drawcountries(ax=ax, linewidth=0.05)

    pcm = plt.pcolormesh(lon_ts, lat_ts, field2plot, norm=norm, cmap=cmap)
#    cbar = plt.colorbar(pcm, norm=norm, cmap=cmap, orientation='vertical', pad=0.05, aspect=15,
#                            shrink=0.7, extend='both')
#    cbar.set_ticks(bounds)
    figname = str(findex).zfill(5)
    plt.savefig(figdir+figname, facecolor="black", dpi=150)
    plt.close()

    print findex

    return findex


def load_wmop_variable(datafile, variable, tindex, m):

#   Load variables and coordinates
    f = cf.read(datafile)
    field2plot = f.select(variable)[0]
    lon = field2plot[0].coord('lon').array
    lat = field2plot[0].coord('lat').array
    field2plot = field2plot.array[tindex]

    lon_ts, lat_ts = m(lon, lat)

    u = f.select('eastward_sea_surface_velocity')[0].array[tindex]
    v = f.select('northward_sea_surface_velocity')[0].array[tindex]
    lon = f.select('eastward_sea_surface_velocity')[0].coord('lon').array
    lat = f.select('eastward_sea_surface_velocity')[0].coord('lat').array

    llon, llat = np.meshgrid(lon, lat)
    lon_uv, lat_uv = m(llon, llat)

    print field2plot.shape

    return field2plot, lon_ts, lat_ts

for year in range(yearmin, yearmax+1):
    for month in range(monthmin, monthmax+1):
        yy = str(year)
        mm = str(month).zfill(2)
        ndays = calendar.monthrange(year,month)[1]
        for day in range(1, ndays):
            dd = str(day).zfill(2)

            datafile = 'http://thredds.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop/'+yy+'/'+mm+'/roms_wmop_'+yy+mm+dd+'.nc'
            #datafile = 'http://thredds.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/model_run_aggregation/wmop/runs/wmop_RUN_'+yy+'-'+mm+'-'+dd+'T00:00:00Z'

            print datafile

            for h in range(0, 8):
                findex += 1
                field2plot, lon_ts, lat_ts = load_wmop_variable(datafile, variable, h, m)
                makeplot(field2plot, lon_ts, lat_ts, vmin, vmax, m, findex, shp_info)