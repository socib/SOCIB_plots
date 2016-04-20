#!/usr/bin/python
# # Filter sea-level data
# Filter the sea-level time series from the stations before Ciutadella.<br>
# Input: list of dates for different rissaga events.

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import datetime
import time
import calendar
import matplotlib.text as text
import matplotlib as mpl
import os
from scipy import signal
from matplotlib import dates

rissagadates = [(2016, 2, 7, 11, 13, 0), (2016, 4, 1, 5, 41, 0), (2015, 9, 19, 4, 0, 0), (2015, 8, 1, 2, 21, 0),
                (2015, 7, 30, 15, 56, 0), (2015, 6, 14, 19, 5, 0), (2015, 6, 12, 4, 50, 0), (2015, 5, 6, 3, 3, 0),
                (2015, 5, 5, 23, 4, 0), (2015, 4, 22, 13, 50, 0)]
ndays = 3  # number of days to look before and after
figdir = "/home/ctroupin/Projects/1-Internal/201530627_Rissaga/figures/Filtered/"

# Making the text bigger and the default color for lines black.
mpl.rcParams.update({'font.size': 20})
mpl.rcParams['lines.color'] = 'k'

def main():
    if rissagadates:
        for events in rissagadates:
            filelist, dayrissaga, daybefore, dayafter = create_filelist(events, ndays)
            pressure, pressure_units, pressure_time, time_units = read_data_list(filelist)
            pressure_time2, pressure2 = select_timeperiod(pressure, pressure_time, daybefore, dayafter)
            pressure_filtered, delay, n1 = filter_timeseries(pressure2)
            plot_filtered_timeseries(pressure2, pressure_units, pressure_time2,
                                     time_units, pressure_filtered, delay, figdir, n1, dayrissaga)
    else:
        print("No event specified in the list")

# http://thredds.priv.socib.es/thredds/dodsC/mooring/sea_level/station_santantoni/L1/2016/dep20150305_station-santantoni_L1_2016-04.nc.html
# http://thredds.priv.socib.es/thredds/dodsC/mooring/sea_level/station_pollensa/L1/2016/dep20110701_station-pollensa_L1_2016-04.nc

def create_filelist(events, ndays):
    # Returns a list of files to process based on the date of an event (tuple)
    # Should be obtained using [Data Discovery](http://appstest.socib.es/DataDiscovery/index.jsp)
    filelist = []
    dayrissaga = datetime.datetime(events[0], events[1], events[2], events[3], events[4], events[5])
    daybefore = dayrissaga - datetime.timedelta(days=ndays)
    dayafter = dayrissaga + datetime.timedelta(days=ndays)
    filelist.append((
        "http://thredds.socib.es/thredds/dodsC/mooring/sea_level/"
        "station_pollensa/L1/%s/dep20110701_station-pollensa_L1_%s-%s.nc"
        % (str(daybefore.year), str(daybefore.year), str(daybefore.month).zfill(2))
    ))
    filelist.append((
        "http://thredds.socib.es/thredds/dodsC/mooring/sea_level/"
        "station_pollensa/L1/%s/dep20110701_station-pollensa_L1_%s-%s.nc"
        % (str(dayafter.year), str(dayafter.year), str(dayafter.month).zfill(2))
    ))
    filelist = list(set(filelist))
    return filelist, dayrissaga, daybefore, dayafter


def read_data_list(filelist):
    # Return pressure, pressure units, time and time units
    # using the file list
    pressure_time, pressure = np.array([]), np.array([])
    pressure_units, time_units = '', ''
    for file2load in filelist:
        print("Working on %s" % (file2load) )
        with netCDF4.Dataset(file2load) as nc:
            time_units = nc.variables['time'].units
            pressure_units = nc.variables['WTR_PRE'].units
            pressure_time = np.append(pressure_time, nc.variables['time'][:])
            pressure = np.append(pressure, nc.variables['WTR_PRE'][:])
    return pressure, pressure_units, pressure_time, time_units


def select_timeperiod(pressure, pressure_time, daybefore, dayafter):
    # Selection of good time period
    tmin = (daybefore - datetime.datetime(1970, 1, 1)).total_seconds()
    tmax = (dayafter - datetime.datetime(1970, 1, 1)).total_seconds()
    goodtime = np.where((pressure_time >= tmin) & (pressure_time <= tmax))[0]
    return pressure_time[goodtime], pressure[goodtime]


def filter_timeseries(pressure):
    # For the filter, it is necessary to set:
    #  * the cutoff frequency,
    #  * the filter length,
    #  * the window applied to the filter
    n1 = 128  # filter length
    windowsname = 'blackman'  # windows name
    sample_rate = 1 / 60.0  # one data per minute
    nyq_rate = sample_rate / 2.0  # Nyquist rate
    frqcut = 1.0 / (36 * 3600.)  # cutoff frequency
    cutoff = frqcut / nyq_rate  # Cutoff relative to Nyquist rate
    taps1 = signal.firwin(n1, cutoff=cutoff, window=(windowsname))
    pressure_filtered = signal.lfilter(taps1, 1.0, pressure)
    delay = 0.5 * (n1 - 1) / sample_rate
    return pressure_filtered, delay, n1


def plot_filtered_timeseries(pressure, pressure_units, pressure_time, time_units,
                             pressure_filtered, delay, figdir, n1, dayrissaga):
    # Very simple plot to show the evolution of the pressure
    hfmt = dates.DateFormatter('%d %B')
    time2plot = netCDF4.num2date(pressure_time, time_units)
    time2plot_filter = netCDF4.num2date(pressure_time - delay, time_units)
    fig = plt.figure(num=None, figsize=(14, 6))
    ax = fig.add_subplot(111)
    plt.plot(time2plot, pressure, 'k', lw=0.5, label='Raw signal')
    plt.plot(time2plot_filter[n1:], pressure_filtered[n1:], 'c', linewidth=2, zorder=2, label='Filtered signal')
    plt.axvline(x=dayrissaga, linewidth=3, color='r', alpha=0.5)
    plt.xlabel('Time')
    plt.ylabel(("Pressure\n (%s)" % (pressure_units)), ha='right', rotation=0)
    plt.legend()
    ax.xaxis.set_major_locator(dates.DayLocator())
    ax.xaxis.set_major_formatter(hfmt)
    plt.grid()
    plt.savefig(os.path.join(figdir, 'SantAntoni_timeseries_' + dayrissaga.strftime('%Y%m%d')))
    plt.close()

    fig = plt.figure(num=None, figsize=(14, 6))
    ax = fig.add_subplot(111)
    plt.plot(time2plot[n1 / 2:-n1 / 2], pressure[n1 / 2:-n1 / 2] - pressure_filtered[n1:], 'k', lw=0.5)
    plt.axvline(x=dayrissaga, linewidth=3, color='r', alpha=0.5)
    plt.xlabel('Time')
    plt.ylabel(("Pressure anomaly\n (%s)" % pressure_units), ha='right', rotation=0)
    ax.set_xlim(time2plot[0], time2plot[-1])
    ax.xaxis.set_major_locator(dates.DayLocator())
    ax.xaxis.set_major_formatter(hfmt)
    fig.autofmt_xdate()
    plt.grid()
    plt.savefig(os.path.join(figdir, 'SantAntoni_anomalies_' + dayrissaga.strftime('%Y%m%d')))
    plt.close()

if __name__ == '__main__':
    main()
