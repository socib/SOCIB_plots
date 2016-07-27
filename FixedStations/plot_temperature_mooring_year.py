#!/usr/bin/python
# coding: utf-8

# # Plot sea water temperature at Buoy "Bahia de Palma"
# The temperature time series is represented for selected months and years, the months are overlaid in order to make easier the comparison.

# In[1]:

import glob
import os
import shutil
import numpy as np
import netCDF4
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
from matplotlib import dates
import datetime, time, calendar
import locale
import matplotlib.font_manager as fm
import matplotlib.image as image
locale.setlocale(locale.LC_ALL, 'en_US.utf8')
mpl.rcParams.update({'font.size': 20})
prop = fm.FontProperties(fname='/home/ctroupin/.fonts/Cube-Regular2.ttf')

# In[2]:
timenow = datetime.datetime.now().strftime('%Y%m%d_%H%M')
figdir = "/home/ctroupin/Pictures/SOCIB"
figdir2 = "/home/ctroupin/public_html/TemperatureTimeSeries"

figname1 = "temp_bahiadepalma_" + timenow + '.png'
figname1latest = "temp_bahiadepalma_latest.png"

figname2 = "temp_canaldeibiza_" + timenow + '.png'
figname2latest = "temp_canaldeibiza_latest.png"

# In[3]:

file_basename = "http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/"
file_list = ['buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-06.nc',
                'buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-07.nc',
                'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-06.nc',
                'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-07.nc',
                'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-06.nc',
                'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-07.nc']
file_list = [file_basename + s for s in file_list]

file_list2 = ['buoy_canaldeibiza-scb_sbe37006/L1/2014/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2014-07.nc',
             'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-06.nc',
             'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-07.nc',
             'buoy_canaldeibiza-scb_sbe37006/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2016-06.nc',
             'buoy_canaldeibiza-scb_sbe37006/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2016-07.nc',
             'buoy_canaldeibiza-scb_sbe37005/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37005_L1_2016-07.nc']
file_list2 = [file_basename + s for s in file_list2]

# # Load data

# In[4]:

def read_time_temp_mooring(filelist):
    buoytemperature = np.array([])
    buoytemperature_QC = np.array([])
    buoytime = np.array([])
    buoydate = np.array([])
    buoyyear = np.array([])

    for datafiles in filelist:
        # print(datafiles)
        with netCDF4.Dataset(datafiles) as nc:
            buoytemperature_QC = np.hstack((buoytemperature_QC, nc.variables['QC_WTR_TEM_SBE37'][:]))
            buoytemperature = np.hstack((buoytemperature, nc.variables['WTR_TEM_SBE37'][:]))
            ttime = nc.variables['time'][:]
            buoytime = np.hstack((buoytime, ttime))
            buoytimeunits = nc.variables['time'].units
            buoyyear = np.hstack((buoyyear,
                                  netCDF4.num2date(ttime[0], buoytimeunits).year * np.ones_like(ttime)))
    buoytemperature = np.ma.masked_where(buoytemperature_QC !=1, buoytemperature)
    buoydate = netCDF4.num2date(buoytime, buoytimeunits)
    return buoytemperature, buoytime, buoydate, buoyyear

def plot_mooring_timeseries(buoyyear, buoydate, buoytemperature, figtitle, figname):
    colorlist = ['c', 'y', 'b']
    hfmt = dates.DateFormatter('%B')
    im = image.imread('/home/ctroupin/Presentations/figures4presentations/logo/logo_socib_square.png')

    fig, ax= plt.subplots(num=None, figsize=(15, 8))

    i=0
    for years in np.unique(buoyyear):
        #print years
        indices = np.where(buoyyear == years)[0]
        #print(indices.min(), indices.max())
        plt.plot(buoydate[indices],
                 buoytemperature[indices],
                 marker = 'o',
                 markeredgecolor = colorlist[i],
                 color = colorlist[i], ms=1, label=int(years))
        i += 1
    ax.xaxis.set_major_locator(dates.MonthLocator())
    ax.xaxis.set_minor_locator(dates.DayLocator())
    ax.xaxis.set_major_formatter(hfmt)
    plt.title(figtitle, fontproperties=prop)
    #fig.autofmt_xdate()
    hl = plt.legend(loc=2, prop=prop)
    plt.grid()
    newax = fig.add_axes([0.8, 0.125, 0.1, 0.1], anchor='NE', zorder=1)
    newax.imshow(im)
    newax.axis('off')
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop)
    plt.savefig(figname, dpi=300)
    #plt.show()
    plt.close()

# Bahia de Palma
figtitle = 'Sea water temperature ($^{\circ}$C)\n at Bahia de Palma buoy'
buoytemperature, buoytime, buoydate, buoyyear = read_time_temp_mooring(file_list)
buoydate = np.array([dd.replace(year=int(buoyyear.min())) for dd in buoydate])
plot_mooring_timeseries(buoyyear, buoydate, buoytemperature, figtitle, os.path.join(figdir, figname1))

# Canal de Ibiza
figtitle = 'Sea water temperature ($^{\circ}$C)\n at Canal de Ibiza buoy'
buoytemperature, buoytime, buoydate, buoyyear = read_time_temp_mooring(file_list2)
buoydate = np.array([dd.replace(year=int(buoyyear.min())) for dd in buoydate])
plot_mooring_timeseries(buoyyear, buoydate, buoytemperature, figtitle, os.path.join(figdir, figname2))




# Copy the figure in public html directory
shutil.copy2(os.path.join(figdir, figname1), os.path.join(figdir2, figname1latest))
shutil.copy2(os.path.join(figdir, figname2), os.path.join(figdir2, figname2latest))
