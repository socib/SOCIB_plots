__author__ = 'ctroupin'

import numpy as np
import cf
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import dates

figdir = "/home/ctroupin/Pictures/DataPlots/"
datafile = "http://thredds.socib.es/thredds/dodsC/drifter/surface_drifter/drifter_svp033-ieo_svp005/L1/2014/dep0001_drifter-svp033_ieo-svp005_L1_2014-04-20.nc"

cmap = plt.cm.hot_r
fmin, fmax = 16., 30.
norm = colors.Normalize(vmin=fmin, vmax=fmax)
bounds = np.arange(fmin, fmax, 2.0)

# Prepare the plot
coordinates = np.array((-72, -10., 10, 35.))
dlon, dlat = 10., 5.
m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
            urcrnrlat=coordinates[3], lat_ts=0.5*(coordinates[2]+coordinates[3]), resolution='i')

# Load the data
f = cf.read(datafile)
temperature = f.select('sea_water_temperature')
lon = temperature[0].coord('lon').array
lat = temperature[0].coord('lat').array
time = temperature[0].coord('time').array
temperatureQC = temperature[0].ancillary_variables[0].array
temperature = temperature[0].array
cf.close_files()

temperature = np.ma.masked_where(temperatureQC != 1, temperature)
temperature = np.ma.masked_where(temperature > 32., temperature)
temperature = np.ma.masked_where(temperature < 18., temperature)

print temperature.min()
with netCDF4.Dataset(datafile) as nc:
    time2plot = netCDF4.num2date(nc.variables['time'][:], nc.variables['time'].units)

temperature = np.ma.masked_where(lon<coordinates[0], temperature)
lon = np.ma.masked_where(lon<coordinates[0], lon)
lat = np.ma.masked_where(lon<coordinates[0], lat)
time = np.ma.masked_where(lon<coordinates[0], time)
lon, lat = m(lon, lat)

fig = plt.figure(figsize=(15, 6))

ax = fig.add_subplot(122)
monthformat = dates.DateFormatter('%B')
monthticks  = dates.MonthLocator(interval=2)  # every month
dayxticks  = dates.DayLocator()  # every month
ax.xaxis.set_major_locator(monthticks)
ax.xaxis.set_major_formatter(monthformat)
#ax.xaxis.set_minor_locator(dayxticks)
plt.scatter(time2plot, temperature, s=5, c=temperature, edgecolor='None', cmap=cmap, norm=norm)
plt.plot(time2plot, temperature, 'k-', lw=0.2)
fig.autofmt_xdate()
plt.ylim(17, 35)
ax.relim()
ax.autoscale_view()

ax = fig.add_subplot(121)
scat = m.scatter(lon, lat, s=10, c=temperature, edgecolor='None', cmap=cmap, norm=norm)
m.drawcoastlines(linewidth=0.2)
m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
                labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)
cbar = plt.colorbar(scat, extend='both', pad=0.05, cmap=cmap, norm=norm, orientation='horizontal')
cbar.set_ticks(bounds)
cbar.set_label('$^{\circ}$C')



plt.savefig(figdir + '20150803_IEO-SVP005_temperature')

plt.show()
plt.close()