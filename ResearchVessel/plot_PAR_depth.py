__author__ = 'ctroupin'


import numpy as np
import cf
#import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
import datetime, time, calendar

datafile = 'http://thredds.socib.es/thredds/dodsC/research_vessel/ctd/socib_rv-scb_sbe9002/L1/2015/dep0011_socib-rv_scb-sbe9002_L1_2015-05-19.nc'

f = cf.read(datafile)
temp = 0.5*(f.select('sea_water_temperature')[0].array+ f.select('sea_water_temperature')[1].array)
psal = 0.5*(f.select('sea_water_practical_salinity')[0].array+ f.select('sea_water_practical_salinity')[1].array)


# Create the figure
fig=plt.figure()
ax = fig.add_subplot(111)

# plt.clabel(cont,inline=True, fmt='%1.1f')
# plt.xlim(smin, smax)
# plt.ylim(tmin, tmax)
plt.plot(psal, temp, 'k.', ms=1.5)
plt.xlabel('Salinity', fontsize=20)
plt.ylabel('Temperature', fontsize=20)
plt.show()
plt.close()