#!/usr/bin/python
__author__ = 'ctroupin'

from numpy import genfromtxt
import matplotlib.pyplot as plt

datafile = "/home/vessel/RTDATA/socib_rv/SCB-MET009/rawArchive/dep0024_2015-07-11_08:00:00/14072015.meteo.proc"
figname = "/home/ctroupin/SOCIB/Facilities/Vessel/SeaBoard/14072015_meteo_proc"

my_data = genfromtxt(datafile, delimiter=',', skip_header=1)

plt.plot(my_data[:,1], 'k', lw=0.5, label='velocidad_real_viento')
plt.plot(my_data[:,2], 'b', lw=0.5, label='velocidad_aparente_viento')
plt.plot(my_data[:,2]- my_data[:,1], 'r', label='difference')
plt.legend(loc=2)
plt.savefig(figname)
plt.show()
plt.close()
