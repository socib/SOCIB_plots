#!/usr/bin/env python
__author__ = 'ctroupin'


import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib import dates
import datetime
import locale

datafile = "http://thredds.socib.es/thredds/dodsC/research_vessel/ctd/"\
    "socib_rv-scb_sbe9002/L1/2015/dep0011_socib-rv_scb-sbe9002_L1_2015-05-19.nc"

figdir = "/home/ctroupin/SOCIB/Facilities/Vessel/PAR/"

with netcdf.Dataset(datafile, 'r', format='NETCDF4') as nc:
    depth = nc.variables['DEPTH'][:]
    SPAR = nc.variables['NET_RAD'][:]
    SPARname = nc.variables['NET_RAD'].getncattr('long_name')
    SPARunits = nc.variables['NET_RAD'].getncattr('units')

nprofiles = depth.shape[0]
print "Number of profiles = " + str(nprofiles)

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(0, nprofiles-1):
    plt.plot(SPAR[i, :], -depth[i, :])

plt.grid()
plt.ylabel("Depth \n (m) ", fontsize=16, rotation=0, horizontalalignment='right')
plt.xlabel(SPARname + "\n" + "(" + SPARunits + ")", fontsize=16)
plt.savefig(figdir + 'SPAR_May2015', dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

plt.close()