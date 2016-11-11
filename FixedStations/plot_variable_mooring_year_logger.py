#!/usr/bin/python
# coding: utf-8

'''
Plot sea water temperature and salinity at "Bahia de Palma" and "Canal de
Ibiza buoys.
The time series are represented for selected months and years,
the months are overlaid in order to make easier the comparison for
consecutive years.
'''

import glob
import os
import shutil
import numpy as np
import netCDF4
import matplotlib as mpl
mpl.use('Agg')      # Necessary for running it with crontab
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
from matplotlib import dates
import datetime, time, calendar
import locale
import matplotlib.font_manager as fm
import matplotlib.image as image
import logging

locale.setlocale(locale.LC_ALL, 'en_US.utf8')
mpl.rcParams.update({'font.size': 20})
# Use the famous SOCIB font: cube!
prop = fm.FontProperties(fname='/home/ctroupin/.fonts/Cube-Regular2.ttf')

def configure_logging():
    logger = logging.getLogger("timeseries_logger")
    logger.setLevel(logging.DEBUG)
    # Format for our loglines
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # Setup console logging
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # Setup file logging as well
    fh = logging.FileHandler('/home/ctroupin/logs/timeseriesplot.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def read_time_temp_mooring(filelist):

    '''
    Read the variables (T and S), time, date and year from a list of files
    '''
    logger = logging.getLogger("timeseries_logger")

    buoytemperature = np.array([])
    buoytemperature_QC = np.array([])
    buoysalinity = np.array([])
    buoysalinity_QC = np.array([])
    buoytime = np.array([])
    buoydate = np.array([])
    buoyyear = np.array([])

    for datafiles in filelist:

        logger.info('Working on {0}'.format(datafiles))
    	try:
            with netCDF4.Dataset(datafiles) as nc:
                buoytemperature_QC = np.hstack((buoytemperature_QC, nc.variables['QC_WTR_TEM_SBE37'][:]))
                buoytemperature = np.hstack((buoytemperature, nc.variables['WTR_TEM_SBE37'][:]))
                buoysalinity_QC = np.hstack((buoysalinity_QC, nc.variables['QC_SALT_SBE37'][:]))
                buoysalinity = np.hstack((buoysalinity, nc.variables['SALT_SBE37'][:]))
                ttime = nc.variables['time'][:]
                buoytime = np.hstack((buoytime, ttime))
                buoytimeunits = nc.variables['time'].units
                buoyyear = np.hstack((buoyyear,
                                      netCDF4.num2date(ttime[0], buoytimeunits).year * np.ones_like(ttime)))
        except RuntimeError:
            logger.error('File {0} does not exist (yet)'.format(datafiles.split('/')[-1]))

    logger.info('Applying mask to the data')
    buoytemperature = np.ma.masked_where(buoytemperature_QC !=1, buoytemperature)
    buoysalinity = np.ma.masked_where(buoysalinity_QC !=1, buoysalinity)
    buoydate = netCDF4.num2date(buoytime, buoytimeunits)
    return buoytemperature, buoysalinity, buoytime, buoydate, buoyyear

def plot_mooring_timeseries(buoyyear, buoydate, buoyvariable, figtitle, figname):

    '''
    Create the figures (time series) for the different variables
    '''

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
                 buoyvariable[indices],
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
    newax = fig.add_axes([0.8, 0.9, 0.1, 0.1], anchor='NE', zorder=1)
    newax.imshow(im)
    newax.axis('off')
    for label in ax.get_xticklabels():
        label.set_fontproperties(prop)

    for label in ax.get_yticklabels():
        label.set_fontproperties(prop)
    plt.savefig(figname, dpi=300)
    #plt.show()
    plt.close()


def main():

    logger = configure_logging()
    logger.info('---Starting new run---')

    # File and directory names
    timenow = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    figdir = "/home/ctroupin/Pictures/SOCIB"
    figdir2 = "/home/ctroupin/public_html/TemperatureTimeSeries"

    figname1T = "temp_bahiadepalma_" + timenow + '.png'
    figname1Tlatest = "temp_bahiadepalma_latest.png"
    figname2T = "temp_canaldeibiza_" + timenow + '.png'
    figname2Tlatest = "temp_canaldeibiza_latest.png"
    figname1S = "psal_bahiadepalma_" + timenow + '.png'
    figname1Slatest = "psal_bahiadepalma_latest.png"
    figname2S = "psal_canaldeibiza_" + timenow + '.png'
    figname2Slatest = "psal_canaldeibiza_latest.png"

    '''
    Generate the file list
    probably a more clever way to do it, but there are several deployments
    for each of the platform and so the names have changed.
    Using new API would be better
    '''

    file_basename = "http://thredds.socib.es/thredds/dodsC/mooring/conductivity_and_temperature_recorder/"
    # file_list = ['buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-06.nc',
    #              'buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-07.nc',
	# 	         'buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-08.nc',
	# 	         'buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-09.nc',
    #              'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-06.nc',
    #              'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-07.nc',
    #              'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-08.nc',
    #              'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-09.nc',
	# 	         'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-06.nc',
    #              'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-07.nc',
    #              'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-08.nc',
	# 	         'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-09.nc']

    file_list = ['buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-10.nc',
                 'buoy_bahiadepalma-scb_sbe37005/L1/2014/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2014-11.nc',
                 'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-10.nc',
                 'buoy_bahiadepalma-scb_sbe37005/L1/2015/dep0002_buoy-bahiadepalma_scb-sbe37005_L1_2015-11.nc',
                 'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-10.nc',
		         'buoy_bahiadepalma-scb_sbe37007/L1/2016/dep0001_buoy-bahiadepalma_scb-sbe37007_L1_2016-11.nc']

    file_list = [file_basename + s for s in file_list]

    # file_list2 = ['buoy_canaldeibiza-scb_sbe37006/L1/2014/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2014-07.nc',
    #               'buoy_canaldeibiza-scb_sbe37006/L1/2014/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2014-08.nc',
	# 	          'buoy_canaldeibiza-scb_sbe37006/L1/2014/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2014-09.nc',
    #               'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-06.nc',
    #               'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-07.nc',
    #               'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-08.nc',
	# 	          'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-09.nc',
    #               'buoy_canaldeibiza-scb_sbe37006/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2016-06.nc',
    #               'buoy_canaldeibiza-scb_sbe37006/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2016-07.nc',
    #               'buoy_canaldeibiza-scb_sbe37005/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37005_L1_2016-08.nc',
	# 	          'buoy_canaldeibiza-scb_sbe37005/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37005_L1_2016-09.nc']

    file_list2 = ['buoy_canaldeibiza-scb_sbe37006/L1/2014/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2014-10.nc',
                  'buoy_canaldeibiza-scb_sbe37006/L1/2014/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2014-11.nc',
                  'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-10.nc',
                  'buoy_canaldeibiza-scb_sbe37006/L1/2015/dep0003_buoy-canaldeibiza_scb-sbe37006_L1_2015-11.nc',
                  'buoy_canaldeibiza-scb_sbe37005/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37005_L1_2016-10.nc',
		          'buoy_canaldeibiza-scb_sbe37005/L1/2016/dep0003_buoy-canaldeibiza_scb-sbe37005_L1_2016-11.nc']

    file_list2 = [file_basename + s for s in file_list2]

    # Bahia de Palma
    mooring = "Bahia de Palma"
    logger.info('Working on %s data' %(mooring))
    figtitleT = 'Sea water temperature ($^{\circ}$C)\n at %s buoy' %(mooring)
    figtitleS = 'Sea water salinity\n at %s buoy' %(mooring)
    logger.debug('Reading the files')
    buoytemperature, buoysalinity, buoytime, buoydate, buoyyear = read_time_temp_mooring(file_list)
    buoydate = np.array([dd.replace(year=int(buoyyear.min())) for dd in buoydate])
    logger.debug('Creating the plots')
    logger.debug('Temperature')
    plot_mooring_timeseries(buoyyear, buoydate, buoytemperature, figtitleT, os.path.join(figdir, figname1T))
    logger.debug('Salinity')
    plot_mooring_timeseries(buoyyear, buoydate, buoysalinity, figtitleS, os.path.join(figdir, figname1S))

    # Canal de Ibiza
    mooring = "Canal de Ibiza"
    logger.info('Working on %s data' %(mooring))
    figtitleT = 'Sea water temperature ($^{\circ}$C)\n at %s buoy' %(mooring)
    figtitleS = 'Sea water salinity\n at %s buoy' %(mooring)
    logger.debug('Reading the files')
    buoytemperature, buoysalinity, buoytime, buoydate, buoyyear = read_time_temp_mooring(file_list2)
    buoydate = np.array([dd.replace(year=int(buoyyear.min())) for dd in buoydate])
    logger.debug('Creating the plots')
    logger.debug('Temperature')
    plot_mooring_timeseries(buoyyear, buoydate, buoytemperature, figtitleT, os.path.join(figdir, figname2T))
    logger.debug('Salinity')
    plot_mooring_timeseries(buoyyear, buoydate, buoysalinity, figtitleS, os.path.join(figdir, figname2S))

    # Copy the figures in public html directory
    logger.info('Making copies of figures into %s' %figdir2)
    shutil.copy2(os.path.join(figdir, figname1T), os.path.join(figdir2, figname1Tlatest))
    shutil.copy2(os.path.join(figdir, figname1S), os.path.join(figdir2, figname1Slatest))
    shutil.copy2(os.path.join(figdir, figname2T), os.path.join(figdir2, figname2Tlatest))
    shutil.copy2(os.path.join(figdir, figname2S), os.path.join(figdir2, figname2Slatest))

if __name__ == "__main__":
    main()
