__author__ = 'ctroupin'

import glob
import netCDF4
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import dates
import datetime
from calendar import monthrange

doplot = 0
plottimeseries = 1
figdir = "/home/ctroupin/DataOceano/MyOcean/figures/mooring/"
mooringdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/history/mooring/2plot/"
mooringlist = sorted(glob.glob(mooringdir+'*.nc'))

def time2decimalday(timearray):
    decimalday = []
    for tt in timearray:
        dayint = tt.timetuple().tm_yday
        daydec = (tt.timetuple().tm_sec + tt.timetuple().tm_min*60. + tt.timetuple().tm_hour*3600.)/86400.
        decimalday.append(dayint+daydec-1)

    decimalday = np.array(decimalday)
    return decimalday

newxticks = []
monthname = []
for m in range(1, 13):
    newxticks.append((datetime.date(2015, m, 1) - datetime.date(2015, 1, 1)).days)
    monthname.append((datetime.date(2015, m, 1).strftime('%B')))

for moorings in mooringlist:
    print moorings

    with netCDF4.Dataset(moorings) as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]
        ttime = nc.variables['TIME']
        time2plot = netCDF4.num2date(ttime[:], ttime.units)
        timevec = ttime[:]

        temp = nc.variables['TEMP'][:]
        tempQC = nc.variables['TEMP_QC'][:]
        temp2 = np.copy(temp)
        temp2 = np.ma.masked_where(tempQC!=1, temp2)

    day2015 = (datetime.datetime(2015,1,1, 0, 0) - datetime.datetime(1950,1,1, 0, 0)).total_seconds()/86400.
    yearindex2015 = np.where(timevec>=day2015)[0]
    yearindex = np.where(timevec<day2015)[0]
    print day2015
    print yearindex

#    fig = plt.figure(figsize=(15, 6))
#    ax = fig.add_subplot(111)
#    plt.plot(time2plot, temp, lw=1)
#    plt.plot(time2plot, temp2, lw=1)

#     plt.title(str(lat.mean()) + 'N - ' + str(lon.mean()) +'E', fontsize=24)
#     fig.autofmt_xdate()
#     ax.set_ylim(12.5, 30.)
#     figname = os.path.basename(moorings).split('.')[0]
# #       plt.savefig(figdir + figname)
#     #plt.show()
#     plt.close()

    day2plot = time2decimalday(time2plot)
    daymask = np.ma.make_mask_none(temp2.shape)
    print day2plot.min()
    print day2plot.max()
    temp3 = np.ma.array(temp2)
    # Compute min, max and mean for each day of the year
    tempmin, tempmax, tempmean, tempstd = [], [], [], []


    for dd in range(0, 366):
        dayindex0 = np.where(np.logical_and((day2plot>=dd), (day2plot<dd+1)))[0]
        dayindex = np.intersect1d(dayindex0, yearindex)
        temp2calc = temp2[dayindex]
        tempstd.append(temp2calc.std())

        # Remove measurements outside 3 standards deviations
        goodvalues = np.where(abs(temp2calc-temp2calc.mean())>=2*temp2calc.std())[0]
        daymask[dayindex[goodvalues]] = True

        temp3 = np.ma.masked_array(temp2, daymask)
        temp2calc = temp3[dayindex]
        tempmin.append(temp2calc.min())
        tempmax.append(temp2calc.max())
        tempmean.append(temp2calc.mean())

    print '**************' + str(len(np.nonzero(daymask)[0]))

    tempmean = np.array(tempmean)
    tempstd = np.array(tempstd)

    print tempmean.shape
    print tempstd.shape

    for dd in range(0, 200):
        dayindex0 = np.where(np.logical_and((day2plot>=dd), (day2plot<dd+1)))[0]
        dayindex = np.intersect1d(dayindex0, yearindex2015)

        temp2calc = temp2[dayindex]

        # Remove measurements outside 3 standards deviations
        goodvalues = np.where(abs(temp2calc-tempmean[dd])>=3*tempstd[dd])[0]

        daymask[dayindex[goodvalues]] = True

        temp3 = np.ma.masked_array(temp2, daymask)





    print temp2.mask
    print temp2.min()
    print temp2.max()
    print '**************' + str(len(np.nonzero(daymask)[0]))
    print np.sum(np.isnan(temp2))
    print 'Min temp = ' + str(np.nanmin(temp2))

    fig2 = plt.figure(figsize=(15, 6))
    ax2 = fig2.add_subplot(111)

    #plt.plot(range(0, 366), tempmin, 'g', lw=2)
    #plt.plot(range(0, 366), tempmax, 'r', lw=2)
    #plt.plot(range(0, 366), tempmean, color='0.25', lw=2, zorder=3)
    #ax2.fill_between(range(0, 366), tempmin, tempmax,facecolor='0.75', zorder=2, alpha=0.5)
    ax2.fill_between(range(0, 366), tempmean-2*tempstd, tempmean+2*tempstd,facecolor='0.25', zorder=2, alpha=0.25)
    #plt.plot(day2plot[yearindex], temp2[yearindex], 'ko', ms=1, markerfacecolor='0.5', lw=0.5)
    # plt.plot(day2plot[yearindex2015], temp2[yearindex2015], 'ro', markerfacecolor='r', ms=2, lw=2)
    plt.plot(day2plot[yearindex2015], temp3[yearindex2015], 'bo', markerfacecolor='b', ms=2, lw=2)

    ax2.set_xlim(0, 365)
    #plt.xticks(newxticks, monthname)
    #fig2.autofmt_xdate()
    plt.ylabel('Sea water \n temperature \n ($^{\circ}$C)', rotation=0, ha='right')
#        plt.savefig(figdir + 'tetetet')
    plt.show()
    plt.close()
