#!/usr/bin/env python

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import dates
import datetime
import locale
import os

locale.setlocale(locale.LC_ALL, 'en_US.utf8')

figdir = '/home/ctroupin/Pictures/'
hfmt = dates.DateFormatter('%d')

month = 1
yearstart, yearend = 2014, 2016

fig = plt.figure(num=None, figsize=(16, 6))
ax = fig.add_subplot(111)

for year in range(yearstart, yearend+1):
    datafile = "http://thredds.socib.es/thredds/dodsC/mooring/weather_station/station_parcbit-scb_met004/L1/" + \
    str(year) + "/dep0002_station-parcbit_scb-met004_L1_" + str(year) +  "-" + str(month).zfill(2) + ".nc"
    with netCDF4.Dataset(datafile) as nc:
        temperature = nc.variables['AIR_TEM'][:]
        temperatureQC = nc.variables['QC_AIR_TEM'][:]
        ttime = nc.variables['time'][:]

    temperature = np.ma.masked_where(temperatureQC != 1, temperature)

    # matplotlib date format object
    dts = map(datetime.datetime.fromtimestamp, ttime - (year-yearstart)*365.*86400)
    fds1 = dates.date2num(dts)

    if year == yearstart:
        plt.xlim(fds1[0]-0.5, fds1[-1]+0.5)

    alpha = 0.65
    if year == yearend:
        alpha = 1
    plt.plot(fds1, temperature, ms=1, lw=1.5, label=str(year), alpha=alpha)


ax.xaxis.set_major_locator(dates.DayLocator(interval=3))
ax.xaxis.set_major_formatter(hfmt)

ax.legend(loc='upper center')
plt.ylabel('$^{\circ}$C', rotation=0, ha='right', fontsize=18)
plt.xlabel('January', fontsize=18)
plt.title('Air temperature at ParcBit')
#fig.autofmt_xdate()
plt.grid()
plt.savefig(figdir + 'temp_parcbit_2016', dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
plt.close()

