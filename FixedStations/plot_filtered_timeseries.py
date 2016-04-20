#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# plot_filtered_timeseries
#
# ctroupin, November 2014
# -------------------------------------------------------------

import os
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib  import colors
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
from scipy import signal

datafile='http://thredds.socib.es/thredds/dodsC/mooring/weather_station/station_parcbit-scb_met004/L1/dep0002_station-parcbit_scb-met004_L1_latest.nc'
with netcdf.Dataset(datafile) as nc:
    AIR_TEM=nc.variables['AIR_TEM'][:]
    time=nc.variables['time'][:]

# Filter AIR_TEM
N=1440
AIR_TEM_filtered=np.convolve(AIR_TEM, np.ones((N,))/N, mode='same')

fig=plt.figure(num=None, figsize=(10,6))
plt.plot(time,AIR_TEM,'k',lw=0.5)
plt.plot(time,AIR_TEM_filtered,'k',lw=2.5)
plt.show()
plt.close()
