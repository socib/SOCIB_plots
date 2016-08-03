__author__ = 'ctroupin'


coriolisdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/profiler-glider/socib_imedea_files/"

datafilesocib_dt = "http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0011_sdeep00_scb-sldeep000_L1_2014-04-07_data_dt.nc"
datafilesocib_rt = "http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0011_sdeep00_scb-sldeep000_L1_2014-04-07_data_rt.nc"
datafilecoriolis = coriolisdir + "GL_201404_PR_GL_68457.nc"

#datafilesocib_dt = "http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0010_sdeep00_scb-sldeep000_L1_2014-02-06_data_dt.nc"
#datafilesocib_rt= "http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0010_sdeep00_scb-sldeep000_L1_2014-02-06_data_rt.nc"
#datafilecoriolis = coriolisdir + "GL_201402_PR_GL_68457.nc"

#datafilesocib_dt = "http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2014/dep0013_ideep00_ime-sldeep000_L1_2014-10-07_data_dt.nc"
#datafilesocib_rt = "http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2014/dep0013_ideep00_ime-sldeep000_L1_2014-10-07_data_rt.nc"
#datafilecoriolis = coriolisdir + "GL_201410_PR_GL_68452.nc"

#datafilesocib = "http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2015/dep0020_sdeep01_scb-sldeep001_L1_2015-01-28_data_rt.nc"
#datafilecoriolis = coriolisdir + "GL_201501_PR_GL_68967.nc"
#datafilecoriolis2 = coriolisdir + "GL_201502_PR_GL_68967.nc"

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import os

figdir = "/home/ctroupin/Projects/201501_InsTAC/GliderData/ComparisonTracks/"
figname = os.path.basename(datafilecoriolis)[:-3]


coordinates = (-.5, 2.5, 38.25, 40.)
m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
            urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
            lat_ts=0.5*(coordinates[2]+coordinates[3]),
            resolution='h')

    #
    # print 'min. longitude = ' + str(lon.min())
    # print 'max. longitude = ' + str(lon.max())
    # print 'min. latitude = ' + str(lat.min())
    # print 'max. latitude = ' + str(lat.max())
    # print ' '
    #
    # lon = np.ma.masked_outside(lon, -360., 360.)
    # lat = np.ma.masked_outside(lat, -90., 90.)
    #
    # print 'min. longitude = ' + str(lon.min())
    # print 'max. longitude = ' + str(lon.max())
    # print 'min. latitude = ' + str(lat.min())
    # print 'max. latitude = ' + str(lat.max())
    # print ' '

fig = plt.figure()
ax = fig.add_subplot(111)
m.ax = ax

with netCDF4.Dataset(datafilesocib_dt) as nc:
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    print 'Fill value = ' + str(nc.variables['longitude']._FillValue)

ax.set_xlim(lon.min(), lon.max())
ax.set_ylim(lat.min(), lat.max())

lon, lat = m(lon, lat)
m.plot(lon, lat, 'r--', label='SOCIB Delayed Mode')

with netCDF4.Dataset(datafilesocib_rt) as nc:
    lon = nc.variables['longitude'][:]
    lat = nc.variables['latitude'][:]
    print 'Fill value = ' + str(nc.variables['longitude']._FillValue)

print 'min. longitude = ' + str(lon.min())
print 'max. longitude = ' + str(lon.max())
print 'min. latitude = ' + str(lat.min())
print 'max. latitude = ' + str(lat.max())

lon[lon > 1000.] = 0.0
lat[lat > 1000.] = 0.0

lon, lat = m(lon, lat)
m.plot(lon, lat, 'go--', ms=1, label='SOCIB Real-time')

with netCDF4.Dataset(datafilecoriolis, 'r') as nc:
     lon = nc.variables['LONGITUDE'][:]
     lat = nc.variables['LATITUDE'][:]

print 'min. longitude = ' + str(lon.min())
print 'max. longitude = ' + str(lon.max())
print 'min. latitude = ' + str(lat.min())
print 'max. latitude = ' + str(lat.max())

lon[lon > 1000.] = 0.0
lat[lat > 1000.] = 0.0

lon, lat = m(lon, lat)
m.plot(lon, lat, 'bo--', ms=1, label='Copernicus')

# with netCDF4.Dataset(datafilecoriolis2, 'r') as nc:
#      lon = nc.variables['LONGITUDE'][:]
#      lat = nc.variables['LATITUDE'][:]
#
# lon, lat = m(lon, lat)
# m.plot(lon, lat, 'bo--', ms=0.5)


m.drawcoastlines()
plt.title(os.path.basename(datafilecoriolis)[:-3], fontsize=20)
plt.legend(loc=2)
plt.savefig(figdir + figname)
# plt.show()
plt.close()
