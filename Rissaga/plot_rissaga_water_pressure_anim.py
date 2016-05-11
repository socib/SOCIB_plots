#!/usr/bin/env python
#
# -------------------------------------------------------------------------

import os
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import datetime, time, calendar

doplot = 1
rissagayear, rissagamonth, rissagaday, rissagahour, rissagadelta = 2016, 4, 1, 8, 2

figdir = ('/home/ctroupin/Projects/1-Internal/201530627_Rissaga/%s%s%s/SeaWaterPressure/' %(str(rissagayear),str(rissagamonth).zfill(2),str(rissagaday).zfill(2)))
figdir
figbasename = 'rissaga_p_'

timerissaga = datetime.datetime(rissagayear, rissagamonth, rissagaday, rissagahour)
trissaga = time.mktime(timerissaga.timetuple())
tmin,tmax = trissaga - rissagadelta * 3600., trissaga + rissagadelta * 3600.

# Generate file name using specified dates
datafile = ('http://thredds.socib.es/thredds/dodsC/mooring/current_profiler/station_ciutadella-ime_awac001/L1/%s/dep0001_station-ciutadella_ime-awac001_L1_%s-%s.nc' %(str(rissagayear),str(rissagayear),str(rissagamonth).zfill(2)))

print("Working on %s" %(datafile))

with netcdf.Dataset(datafile,'r+', format='NETCDF4') as nc:
    
    pressure_time=nc.variables['time'][:]
    goodtime=np.where(np.logical_and((pressure_time<=tmax),(pressure_time>=tmin)))[0]
    pressure=nc.variables['WTR_PRE'][goodtime]
    pressure_time=pressure_time[goodtime]
    

ii=1
if doplot == 1:
    for tt,press in zip(pressure_time,pressure):
        fig=plt.figure(num=None, figsize=(10, 6))
        #fig.patch.set_alpha(0.)
        plt.plot(pressure_time,pressure,'k')
        plt.plot(tt,press,'co',ms=10)
        #plt.title(time.strftime("%H:%M:%S", time.gmtime(tt)))
        plt.axis('off')
        # plt.show()
        print os.path.join(figdir, figbasename + str(ii).zfill(4))
        plt.savefig(os.path.join(figdir, figbasename + str(ii).zfill(4)), dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)
        ii+=1
        plt.close()
