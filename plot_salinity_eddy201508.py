__author__ = 'ctroupin'

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import cf
from matplotlib import colors

figdir = "/home/ctroupin/Pictures/SOCIB/"
modelfile = "http://thredds.priv.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop/2015/08/roms_wmop_20150814.nc"
drifterfile = "http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp023-scb_svp018/L1/2015/dep0001_drifter-svp023_scb-svp018_L1_2015-07-13.nc"

yearmin, yearmax = 2015, 2015
monthmin, monthmax = 8, 8
vmin, vmax, dvar = 36.5, 38.5, 0.5

variable = 'sea_surface_salinity'
cmap = plt.cm.RdYlBu_r

norm = colors.Normalize(vmin=vmin, vmax=vmax)
normvel = colors.Normalize(vmin=0.0, vmax=0.5)
bounds = np.arange(vmin, vmax + .0001, dvar)


coordinates = np.array((-1, 5.001, 36.9, 41.))
dlon, dlat = 1.0, 1.0

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution='h')

findex = 0

landfile = "/data_local/Bathymetry/ne_10m_land"
shp_info = m.readshapefile(landfile, 'scalerank', drawbounds=True)

f = cf.read(drifterfile)
lon_d = f.select('sea_water_temperature')[0].coord('lon').array
lat_d = f.select('sea_water_temperature')[0].coord('lat').array
f.close()

lon_d = np.ma.masked_outside(lon_d, coordinates[0], coordinates[1], copy='True')
print lon_d.min()
print lon_d.max()
lon_d, lat_d = m(lon_d, lat_d)

def load_wmop_variable(datafile, variable, tindex, m):

#   Load variables and coordinates
    f = cf.read(datafile)
    field2plot = f.select(variable)[0]
    lon = field2plot[0].coord('lon').array
    lat = field2plot[0].coord('lat').array
    field2plot = field2plot.array[tindex]

    lon, lat = np.meshgrid(lon, lat)
    lon_ts, lat_ts = m(lon, lat)

    u = f.select('eastward_sea_surface_velocity')[0].array[tindex]
    v = f.select('northward_sea_surface_velocity')[0].array[tindex]
    lon = f.select('eastward_sea_surface_velocity')[0].coord('lon').array
    lat = f.select('eastward_sea_surface_velocity')[0].coord('lat').array

    llon, llat = np.meshgrid(lon, lat)
    lon_uv, lat_uv = m(llon, llat)

    print field2plot.shape

    return field2plot, lon_ts, lat_ts, lon_uv, lat_uv, u, v

field2plot, lon_ts, lat_ts, lon_uv, lat_uv, u, v = load_wmop_variable(modelfile, 'sea_surface_salinity', 12, m)


fig = plt.figure()
ax = fig.add_subplot(111)
ax = plt.gca()
ax.cla()


m.drawcoastlines(ax=ax, linewidth=0.05)
m.drawcountries(ax=ax, linewidth=0.05)

m.fillcontinents(color = 'gray')
m.drawmapboundary(linewidth=0.5, color='k')
m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
                labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)


m.drawmapscale(-0.25, 41, -0.75, 41, 50, barstyle='simple', units='km', fontsize=12, zorder=3)

m.plot(lon_d, lat_d, 'ko', ms=0.5, alpha=0.5)

pcm = m.pcolormesh(lon_ts, lat_ts, field2plot, norm=norm, cmap=cmap)
cbar = plt.colorbar(pcm, norm=norm, cmap=cmap, orientation='vertical', pad=0.05, aspect=15,
                        shrink=0.9, extend='both')
cbar.set_ticks(bounds)

speed = np.sqrt(u.data * u.data + v.data * v.data)
speed = np.ma.masked_greater(speed, 1.5)
plt.streamplot(lon_uv, lat_uv, u, v, color=speed,
                arrowsize=1.5, norm=normvel,cmap=plt.cm.gray_r,
                density=15, linewidth=.5, zorder=4)

plt.title('WMOP Salinity and drifter SVP018\n 2015-08-14')
figname = str(findex).zfill(5)
plt.savefig(figdir+'salinity_drifter20150825_V3', dpi=300)
# plt.show()
plt.close()
