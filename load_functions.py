import os
import numpy as np
import netCDF4 
import cf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png

def read_profiler_data(datafile):
    with netCDF4.Dataset(datafile) as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        time = nc.variables['time'][:]
        depth = nc.variables['DEPTH'][:]
        temp = nc.variables['WTR_TEM'][:]
        psal = nc.variables['SALT_ADJUSTED'][:]
    return lon, lat, depth, time, temp, psal

def load_glider_position(filename):
    with netCDF4.Dataset(filename) as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
    return lon, lat

def load_glider_TS(filename):
    with netCDF4.Dataset(filename) as nc:
	psal = nc.variables['salinity'][:]
	temp = nc.variables['temperature'][:]
    return temp, psal

def load_turtle_coord(filename):
    with netCDF4.Dataset(filename) as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        lonQC = nc.variables['QC_LON'][:]
        latQC = nc.variables['QC_LAT'][:]
	return lon[lonQC == 1], lat[latQC == 1]

def load_profiler_TS(filename):
    f = cf.read(filename)
    temp = f.select('sea_water_temperature')[1]
    psal = f.select('sea_water_salinity')[1]
    f.close()
    return temp.array, psal.array

def load_salinity_L4_SMOS(filename):
    with netCDF4.Dataset(filename) as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        psal = nc.variables['l4_sss'][:].squeeze()
    return lon, lat, psal



def load_profiler_data(filename):
    f = cf.read(filename)
    temp = f.select('sea_water_temperature')[1]
    psal = f.select('sea_water_salinity')[1]
    depth = temp.coord('depth')
    time = temp.coord('time')
    lon = temp.coord('longitude')
    lat = temp.coord('latitude')
    f.close()
    return lon, lat, depth, time, temp.array, psal.array

def load_sst_modis(filename): 
    
    if 'SST4' in os.path.basename(filename):
        var2load = 'sst4'
        varqc2load = 'qual_sst4'
    else:
        var2load = 'sst'
        varqc2load = 'qual_sst'

    with netCDF4.Dataset(filename) as nc:
        lon = nc.groups['navigation_data'].variables['longitude'][:]
        lat = nc.groups['navigation_data'].variables['latitude'][:]
        sst = nc.groups['geophysical_data'].variables[var2load][:]
        qualsst = nc.groups['geophysical_data'].variables[varqc2load][:]

        # apply QC filter
        sst = np.ma.masked_where(qualsst > 1, sst)
    return lon, lat, sst

def load_altimetry_aviso_uv(altimetryfile, coordinates):
    with netCDF4.Dataset(altimetryfile) as nc:
        lon = nc.variables['lon'][:] - 360.
        lat = nc.variables['lat'][:]
        u = np.squeeze(nc.variables['u'][:])
        v = np.squeeze(nc.variables['v'][:])
        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        u = u[goodlat, :]
        u = u[:, goodlon]
        v = v[goodlat, :]
        v = v[:, goodlon]
    return lon, lat, u, v

def load_altimetry_aviso_adt(altimetryfile, coordinates):
    with netCDF4.Dataset(altimetryfile) as nc:
        lon = nc.variables['lon'][:] - 360.
        lat = nc.variables['lat'][:]
        adt = nc.variables['adt'][:].squeeze()
	# subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        adt = adt[goodlat, :]
        adt = adt[:, goodlon]
    return lon, lat, adt
