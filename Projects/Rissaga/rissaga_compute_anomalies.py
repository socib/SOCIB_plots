#!/usr/bin/env python
#
# -------------------------------------------------------------------------

import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib import rcParams
import datetime
import time
from scipy import signal

figdir = '/home/ctroupin/Projects/20153627_Rissaga/figures/'
figbasename = 'rissaga_p_'

rcParams['axes.color_cycle'] = ["c", "b", "r", "g"]

# Define time limits for the plot
timeinit, timeend = datetime.datetime(2015, 6, 1, 00, 00), datetime.datetime(2015, 6, 30, 1, 0)
tmin, tmax = time.mktime(timeinit.timetuple()), time.mktime(timeend.timetuple())

# ciutadellafile='http://thredds.socib.es/thredds/dodsC/mooring/current_profiler/station_ciutadella-ime_awac001/L1/dep0001_station-ciutadella_ime-awac001_L1_latest.nc'
ciutadellafile = "http://thredds.socib.es/thredds/dodsC/mooring/current_profiler/station_ciutadella-ime_awac001/L1/2015/dep0001_station-ciutadella_ime-awac001_L1_2015-06.nc"

with netcdf.Dataset(ciutadellafile, 'r+', format='NETCDF4') as nc:
    pressure_time = nc.variables['time'][:]
    goodtime = np.where(np.logical_and((pressure_time <= tmax), (pressure_time >= tmin)))[0]
    pressure = nc.variables['WTR_PRE'][goodtime]
    pressure_time = pressure_time[goodtime]

# matplotlib date format object
hfmt = dates.DateFormatter('%d %B')

# Design the filter
N1 = 128  # filter length
windowsname = 'blackman'  # windows name

sample_rate = 1 / 60.0  # one data per minute
nyq_rate = sample_rate / 2.0  # Nyquist rate
frqcut = 1.0 / (12 * 3600.)  # 12 hours
cutoff = frqcut / nyq_rate  # Cutoff relative to Nyquist rate

# Build and apply filter
taps1 = signal.firwin(N1, cutoff=cutoff, window=windowsname)
pressure_filtered = signal.lfilter(taps1, 1.0, pressure)
delay = 0.5 * (N1 - 1) / sample_rate

dts = map(datetime.datetime.fromtimestamp, pressure_time)
dts_f = map(datetime.datetime.fromtimestamp, pressure_time - delay)
fds = dates.date2num(dts)  # converted
fds_f = dates.date2num(dts_f)

# Compute pressure anomaly
pressure_anomaly = pressure[:-N1 / 2] - pressure_filtered[N1 / 2:]

# Identify period of Rissaga
rissagatime = np.where(abs(pressure_anomaly[N1 / 2:]) > 0.6)[0]
dtime = 6. * 60.

fig = plt.figure(num=None, figsize=(12, 6))
ax = fig.add_subplot(211)

plt.title('Pressure at 6 meters', rotation=0, fontsize=24)

ax.xaxis.set_major_formatter(hfmt)
plt.plot(fds[:], pressure[:], 'k', lw=0.5)
plt.plot(fds_f, pressure_filtered, 'r', lw=2)
ax.xaxis.set_major_locator(dates.DayLocator())
ax.xaxis.set_major_formatter(hfmt)
ax.set_xlim(fds[N1 / 2], fds[-1])
ax.set_ylim(5.5, 7.25)
for rtime in rissagatime:
    ax.axvspan(fds[rtime - dtime], fds[rtime + dtime], alpha=0.05, color='b')

fig.autofmt_xdate()

ax = fig.add_subplot(212)
plt.title('Pressure anomaly at 6 meters', rotation=0, fontsize=24)
ax.xaxis.set_major_locator(dates.DayLocator())
ax.xaxis.set_major_formatter(hfmt)
ax.set_xlim(fds[N1 / 2], fds[-1])
ax.set_ylim(0.0, 1.0)
for rtime in rissagatime:
    ax.axvspan(fds_f[rtime - dtime], fds_f[rtime + dtime], alpha=0.05, color='b')
plt.plot(fds_f[N1 / 2:], abs(pressure_anomaly), 'k', lw=0.5)
fig.autofmt_xdate()

plt.show()
plt.close()
