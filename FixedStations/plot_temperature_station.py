
# coding: utf-8

# # Plot sea water temperature at Buoy "Andratx"

import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib import dates
import matplotlib as mpl
import calendar
import locale
import datetime
locale.setlocale(locale.LC_ALL, 'en_US.utf8')

yearlist = (2012, 2013, 2014, 2015)
monthlist = (6, 7)

figdir = "/home/ctroupin/Pictures/"

# ## Files to read
fileprefix = 'http://thredds.socib.es/thredds/dodsC/mooring/sea_level/station_andratx/L1/'
stationname = 'Andratx'


fig=plt.figure(num=None, figsize=(12, 6))
ax = fig.add_subplot(111)
ax.set_color_cycle(['k', 'b', 'r', 'g'])
hfmt = dates.DateFormatter('%B')

# Loop on years
for years in yearlist:
    yy = str(years)
    temperature2plot = []
    time2plot = []

    deltatime = calendar.timegm((years, 1, 1, 0, 0, 0)) - calendar.timegm((yearlist[0], 1, 1, 0, 0, 0))


    for months in monthlist:
        mm = str(months).zfill(2)
        file2load = fileprefix + yy + "/dep20110602_station-andratx_L1_" + yy + "-" + mm +".nc"
        print file2load

        with netcdf.Dataset(file2load, 'r', format='NETCDF4') as nc:
            temperature2plot = np.append(temperature2plot, nc.variables['WTR_TEM'][:])
            time2plot = np.append(time2plot, nc.variables['time'][:])
            temperature2plot = np.ma.masked_greater(temperature2plot, 40.)


    dts = map(datetime.datetime.fromtimestamp, time2plot - deltatime)
    fds = dates.date2num(dts)

    plt.plot(fds, temperature2plot, '-', ms=1, label=yy)

ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_minor_locator(dates.DayLocator())
ax.xaxis.set_major_formatter(hfmt)
ax.legend(loc=2)
#ax.set_ylim(16., 30.)
plt.title('Sea water temperature ($^{\circ}$C)\n at '+ stationname)
fig.autofmt_xdate()
plt.grid()
plt.savefig(figdir + stationname + str(yearlist[0]) + '_' + str(yearlist[-1]))
plt.show()
plt.close()


