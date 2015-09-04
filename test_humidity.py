__author__ = 'ctroupin'

from meteo.humidity import *
from meteo.air import *
import numpy as np
import cf
import netCDF4
import matplotlib.pyplot as plt

NN = 24.*60.
figdir = "/home/ctroupin/SOCIB/Facilities/FixedStation/figure/"
datafile = "http://thredds.socib.es/thredds/dodsC/mooring/weather_station/station_parcbit-scb_met004/L1/dep0002_station-parcbit_scb-met004_L1_latest.nc"

f = cf.read(datafile)

pressure = f.select('air_pressure')[0].array
temperature = f.select('air_temperature')[0].array
relativehum = f.select('relative_humidity')[0].array

with netCDF4.Dataset(datafile) as nc:
    time = nc.variables['time']
    time2plot = netCDF4.num2date(time[:], time.units)

TEMP_K_C_OFFSET = 273.15

def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]


def C2K(temp_c):
    return temp_c + TEMP_K_C_OFFSET


def calc_dew_c(rel_humid, temp_c):
    l = np.log(rel_humid / 100.0)
    m = 17.27 * temp_c
    n = C2K(temp_c)
    b = (l + (m / n)) / 17.27
    dew_c = (TEMP_K_C_OFFSET * b) / (1 - b)
    return dew_c


def h2o_ppm(rel_humid, temp_c):
    CA = 8.1332
    CB = 1762.39
    CC = 235.66
    dew_c = calc_dew_c(rel_humid, temp_c)
    curr_temp_part_pressure = CA - (CB / (temp_c + CC))
    amb_part_press_mm_hg = 10 ** curr_temp_part_pressure
    dew_temp_part_pressure = CA - (CB / (dew_c + CC))
    dew_temp_part_press_mm_hg = 10 ** dew_temp_part_pressure
    water_ppm = (dew_temp_part_press_mm_hg / 760.0) * 1000000
    return water_ppm


def h2o_grams_per_cubic_meter(rel_humid, temp_c):
    water_ppm = h2o_ppm(rel_humid, temp_c)
    water_gram_per_cubic_meter = water_ppm * 0.001 * 18.0 / 22.4
    return water_gram_per_cubic_meter

absolutehum1 = h2o_grams_per_cubic_meter(relativehum, temperature)
absolutehum2 = rh2ah(relativehum, pressure[0,:]*100., temperature+ TEMP_K_C_OFFSET)

print absolutehum2.shape
print pressure.shape
print temperature.shape

fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(311)
plt.plot(time2plot[-30*24*60:], temperature[-30*24*60:], 'k', lw=0.5)
#plt.plot(time2plot[-30*24*60:], runningMeanFast(temperature, NN)[-30*24*60:], 'k', lw=2)
plt.ylabel('Temperature\n ($^{\circ}$C)', rotation=0, ha='right')
ax = fig.add_subplot(312)
plt.plot(time2plot[-30*24*60:], relativehum[-30*24*60:], 'k', lw=0.5)
#plt.plot(time2plot[-30*24*60:], runningMeanFast(relativehum, NN)[-30*24*60:], 'k', lw=2)
plt.ylabel('Relative\n humidity\n (%)', rotation=0, ha='right')
ax = fig.add_subplot(313)
plt.plot(time2plot[-30*24*60:], absolutehum1[-30*24*60:], 'k', lw=0.5)
#plt.plot(time2plot[-30*24*60:], runningMeanFast(absolutehum1, NN)[-30*24*60:], 'k', lw=2)
#plt.plot(time2plot, absolutehum2*10., 'r')
plt.ylabel('Absolute\n humidity\n (g/m$^3$)', rotation=0, ha='right')
fig.autofmt_xdate()

plt.savefig(figdir+'relative_humidity')
plt.show()
plt.close()