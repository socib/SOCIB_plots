#!/usr/bin/python

import os
import numpy as np
import netCDF4 as netcdf
import cf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png



def load_glider_position(filename):
    with netcdf.Dataset(filename) as nc:
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

    with netcdf.Dataset(filename) as nc:
        lon = nc.groups['navigation_data'].variables['longitude'][:]
        lat = nc.groups['navigation_data'].variables['latitude'][:]
        sst = nc.groups['geophysical_data'].variables[var2load][:]
        qualsst = nc.groups['geophysical_data'].variables[varqc2load][:]

        # apply QC filter
        sst = np.ma.masked_where(qualsst > 1, sst)
    return lon, lat, sst
