#!/usr/bin/python

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



def load_glider_position(filename):
    with netcdf.netCDF4(filename) as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
    return lon, lat

def load_sst_modis(filename): 
    
    if 'SST4' in os.path.basename(filename):
        var2load = 'sst4'
        varqc2load = 'qual_sst4'
    else:
        var2load = 'sst'
        varqc2load = 'qual_sst'

    with netcdf.netCDF4(filename) as nc:
        lon = nc.groups['navigation_data'].variables['longitude'][:]
        lat = nc.groups['navigation_data'].variables['latitude'][:]
        sst = nc.groups['geophysical_data'].variables[var2load][:]
        qualsst = nc.groups['geophysical_data'].variables[varqc2load][:]

        # apply QC filter
        sst = np.ma.masked_where(qualsst > 1, sst)
    return lon, lat, sst

def load_altimetry_aviso(altimetryfile, coordinates):
    with netcdf.netCDF4(altimetryfile) as nc:
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

