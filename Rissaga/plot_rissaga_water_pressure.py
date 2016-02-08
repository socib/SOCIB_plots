#!/usr/bin/env python
#
# -------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import datetime, time, calendar
import matplotlib.text as text

figdir='/home/ctroupin/Projects/201530627_Rissaga/20160207/SeaWaterPressure/'
figbasename='rissaga_p_'

# Define time limits for the plot
timeinit,timeend= datetime.datetime(2016, 2, 7, 0, 00),datetime.datetime(2016, 2, 8, 0, 0)
tmin,tmax=time.mktime(timeinit.timetuple()),time.mktime(timeend.timetuple())

ciutadellafile='http://thredds.socib.es/thredds/dodsC/mooring/current_profiler/station_ciutadella-ime_awac001/L1/dep0001_station-ciutadella_ime-awac001_L1_latest.nc'


with netcdf.Dataset(datafile,'r+', format='NETCDF4') as nc:
    
    pressure_time=nc.variables['time'][:]
    goodtime=np.where(np.logical_and((pressure_time<=tmax),(pressure_time>=tmin)))[0]
    pressure=nc.variables['WTR_PRE'][goodtime]
    pressure_time=pressure_time[goodtime]
    

ii=1
for tt,press in zip(pressure_time,pressure):
    fig=plt.figure(num=None, figsize=(10, 6))
    #fig.patch.set_alpha(0.)
    plt.plot(pressure_time,pressure,'k')
    plt.plot(tt,press,'co',ms=10)
    #plt.title(time.strftime("%H:%M:%S", time.gmtime(tt)))
    plt.axis('off')
    
    plt.savefig(figdir+figbasename+str(ii).zfill(4), dpi=300, facecolor='w', edgecolor='w',
         transparent=False, bbox_inches='tight', pad_inches=0.1)
    ii+=1
    plt.close()
