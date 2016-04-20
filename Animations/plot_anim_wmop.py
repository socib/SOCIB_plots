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

coordinates = np.array((-6.75, 9.001, 34.75, 45.))
dlon, dlat = 3.0, 2.0

m = Basemap(llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
            urcrnrlon=coordinates[1], urcrnrlat=coordinates[3], resolution='h')

findex = 0

landfile = "/data_local/Bathymetry/ne_10m_land"
shp_info = m.readshapefile(landfile, 'scalerank', drawbounds=True)

def makeplot(field2plot, lon_ts, lat_ts, vmin, vmax, m, findex, shp_info):

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    bounds = np.arange(vmin, vmax + .0001, dvar)

    fig=plt.figure(facecolor='black')
    # ax = fig.add_subplot(111)
    ax = plt.gca()
    ax.cla()

    m.drawcoastlines(ax=ax, linewidth=0.05)
    m.drawcountries(ax=ax, linewidth=0.05)

    m.fillcontinents(color = 'black')
    m.drawmapboundary()
#    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
#                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
#    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
#                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)


# u2plot = u[tt, ::NN, ::NN].squeeze()
# v2plot = v[tt, ::NN, ::NN].squeeze()
# speed = np.sqrt(u2plot*u2plot + v2plot*v2plot)
# speed = np.ma.masked_greater(speed, 1.5)

#    pcm = plt.pcolormesh(lon_ts, lat_ts, field2plot, norm=norm, cmap=cmap)
#    cbar = plt.colorbar(pcm, norm=norm, cmap=cmap, orientation='vertical', pad=0.05, aspect=15,
#                            shrink=0.7, extend='both')
#    cbar.set_ticks(bounds)
    figname = str(findex).zfill(5)
    plt.savefig(figdir+figname, facecolor=fig.get_facecolor(), dpi=150)
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