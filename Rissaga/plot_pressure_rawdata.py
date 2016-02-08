#!/usr/bin/env python

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
from matplotlib import dates
import datetime, time, calendar
from scipy import signal
import matplotlib.text as text
import locale
locale.setlocale(locale.LC_ALL, 'en_US.utf8')
import itertools


datadir = '/home/ctroupin/Projects/20153627_Rissaga/rawdata/2015/05/'
figdir = '/home/ctroupin/Projects/20153627_Rissaga/figures/'
filelist = sorted(glob.glob(datadir+'eth*.ail'))

file0 = filelist[0]

pressure=[0]

for files in filelist:
    with open(files) as f:
        for line in itertools.islice(f, 0, None, 6):
            # Do something with the line
            
            pressure.append(line.split(None, 19)[13])

print len(pressure)


fig=plt.figure()
ax = fig.add_subplot(111)
plt.plot(pressure, lw=0.2)
ax.set_ylim(6.0, 6.8)
plt.xlabel('Time')
plt.ylabel('Pressure at 6 meters')

plt.savefig(figdir + 'pressure_raw_data_201505', dpi=300, facecolor='w', edgecolor='w',
                   transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()
