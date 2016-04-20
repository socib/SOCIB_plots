#!/usr/bin/env python

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
from matplotlib import dates
import datetime, time, calendar
from scipy import signal
import matplotlib.text as text
from matplotlib.path import Path
import matplotlib.patches as patches
import pysocibclient
import locale
locale.setlocale(locale.LC_ALL, 'en_US.utf8')

figdir = '/home/ctroupin/Pictures/'

BoyaPalma201306 = 'http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/buoy_bahiadepalma-scb_sbe37005/L1/2013/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2013-06.nc'
BoyaPalma201307 = 'http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/buoy_bahiadepalma-scb_sbe37005/L1/2013/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2013-07.nc'
BoyaPalma201406 = 'http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-06.nc'
BoyaPalma201407 = 'http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-07.nc'
BoyaPalma201506 = 'http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-06.nc'
BoyaPalma201507 = 'http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-07.nc'

hfmt = dates.DateFormatter('%B')


# 2013
with netcdf.Dataset(BoyaPalma201306,'r', format='NETCDF4') as nc:
        temperature = nc.variables['WTR_TEM_SBE37'][:]
        ttime = nc.variables['time'][:]
with netcdf.Dataset(BoyaPalma201307,'r', format='NETCDF4') as nc:
        temperature = np.append(temperature, nc.variables['WTR_TEM_SBE37'][:])
        ttime = np.append(ttime, nc.variables['time'][:])

# Filter bad values
temperature = np.ma.masked_greater(temperature, 28.)

# matplotlib date format object
dts = map(datetime.datetime.fromtimestamp, ttime + 86400.*365.) 
fds1 = dates.date2num(dts) # converted

fig=plt.figure(num=None, figsize=(12, 6))
ax = fig.add_subplot(111)
#plt.plot(fds1, temperature, 'bo-', ms=1, label='2013')


# 2014
with netcdf.Dataset(BoyaPalma201406,'r', format='NETCDF4') as nc:
        temperature = nc.variables['WTR_TEM_SBE37'][:]
        ttime = nc.variables['time'][:]
with netcdf.Dataset(BoyaPalma201407,'r', format='NETCDF4') as nc:
        temperature = np.append(temperature, nc.variables['WTR_TEM_SBE37'][:])
        ttime = np.append(ttime, nc.variables['time'][:])

# Filter bad values
temperature = np.ma.masked_greater(temperature, 28.)

# matplotlib date format object
hfmt = dates.DateFormatter('%B')
dts = map(datetime.datetime.fromtimestamp, ttime ) 
fds = dates.date2num(dts) # converted


plt.plot(fds, temperature, 'ko-', ms=1, label='2014')


# 2015
with netcdf.Dataset(BoyaPalma201506,'r', format='NETCDF4') as nc:
        temperature = nc.variables['WTR_TEM_SBE37'][:]
        ttime = nc.variables['time'][:]
with netcdf.Dataset(BoyaPalma201507,'r', format='NETCDF4') as nc:
        temperature = np.append(temperature, nc.variables['WTR_TEM_SBE37'][:])
        ttime = np.append(ttime, nc.variables['time'][:])

# Filter bad values
temperature = np.ma.masked_greater(temperature, 28.)

dts = map(datetime.datetime.fromtimestamp, ttime-86400.*365.) 
fds2 = dates.date2num(dts) # converted

plt.plot(fds2, temperature, 'ro-', ms=1, markeredgecolor='r', label='2015')

ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_minor_locator(dates.DayLocator())
ax.xaxis.set_major_formatter(hfmt)
ax.set_xlim(fds2.min()-0.5,fds.max()+1)
ax.set_ylim(20.,28.)
ax.legend(loc=2)
plt.title('Sea water temperature ($^{\circ}$C)\n at Bahia de Palma buoy')
fig.autofmt_xdate()
plt.grid()
plt.savefig(figdir + 'temp_bahiadepalma_20150810', dpi=300, facecolor='w', edgecolor='w',
                       transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
