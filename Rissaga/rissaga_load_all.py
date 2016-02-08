#!/usr/bin/env python

import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib import dates
import datetime
from scipy import signal
import pysocibclient
import locale

locale.setlocale(locale.LC_ALL, 'en_US.utf8')

opendaplist = './Ciutadella.txt'
instrument_type = "current_profiler"
instrument_name = "Ciutadella"

figdir = '/home/ctroupin/Projects/20153627_Rissaga/figures/'
figname1 = 'rissage_timeseries_2.eps'

dtime = 30 * 60.


def create_deployement_list(instrument_type, instrument_name):
    socib_api = pysocibclient.ApiClient()
    alist = socib_api.list_platforms(instrument_type=instrument_type)

    mylist = []
    for deployment in alist:
        if instrument_name in deployment.name:
            for productlist in deployment.product_list:
                if productlist.processing_level == "L1":
                    mylist.append(productlist.opendap)
    return mylist

##if os.path.isfile(opendaplist):
##    # The file exists so we can read it
##    print 'File exists'
##    mylist = []
##    with open(opendaplist, 'r') as f:
##        mylist = f.readlines()        
##else:
mylist = create_deployement_list(instrument_type, instrument_name)
with open(opendaplist, 'wb') as f:
    for item in mylist:
        f.write("%s\n" % item)

# Load all the datapres
pressure_time = []
pressure = []
pressure_QC = []

for ffiles in mylist:
    print ffiles
    with netcdf.Dataset(ffiles, 'r', format='NETCDF4') as nc:
        pressure_time = np.append(pressure_time, nc.variables['time'][:])
        pressure = np.append(pressure, nc.variables['WTR_PRE'][:])
        pressure_QC = np.append(pressure_QC, nc.variables['QC_WTR_PRE'][:])
##        pressure_time=pressure_time[goodtime]

pressure2filter = np.copy(pressure)
pressure[pressure_QC > 1] = np.nan
pressure2filter[pressure_QC > 1] = 0.0
pressure2filter[np.isnan(pressure2filter)] = 0.0

# matplotlib date format object
hfmt = dates.DateFormatter('%B')
dts = map(datetime.datetime.fromtimestamp, pressure_time)
fds = dates.date2num(dts)  # converted

# Design the filter
N1 = 128  # filter length
windowsname = 'blackman'  # windows name

sample_rate = 1 / 60.0  # one data per minute
nyq_rate = sample_rate / 2.0  # Nyquist rate
frqcut = 1.0 / (12 * 3600.)  # 12 hours
cutoff = frqcut / nyq_rate  # Cutoff relative to Nyquist rate

# Build and apply filter
taps1 = signal.firwin(N1, cutoff=cutoff, window=(windowsname))
pressure_filtered = signal.lfilter(taps1, 1.0, pressure2filter)
pressure_filtered[pressure_QC > 1] = np.nan
pressure_filtered[pressure_filtered < 6.] = np.nan

delay = 0.5 * (N1 - 1) / sample_rate
dts2 = map(datetime.datetime.fromtimestamp, pressure_time - delay)
fds2 = dates.date2num(dts2)  # converted

# Identify period of Rissaga
rissagatime = np.where(pressure > 6.8)[0]

# Make the plot
fig = plt.figure(num=None, figsize=(12, 6))
ax = fig.add_subplot(111)
ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_minor_locator(dates.DayLocator())
ax.xaxis.set_major_formatter(hfmt)

# Make the plot
plt.plot(fds, pressure, lw=0.5)
# plt.plot(fds2[:], pressure_filtered[:], color='b', lw=1)

for rtime in rissagatime:
    ax.axvspan(fds[rtime - dtime], fds[rtime + dtime], alpha=0.05, color='r')

plt.title('Pressure at 6 meters', rotation=0, fontsize=24)
plt.ylabel('dbar       \n   ', rotation=0)
ax.set_ylim(5.8, 7.0)
fig.autofmt_xdate()
plt.grid()
plt.savefig(figdir + figname1, dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.show()
plt.close()
