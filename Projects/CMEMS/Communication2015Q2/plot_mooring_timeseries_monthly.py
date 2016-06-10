__author__ = 'ctroupin'

import glob
import os
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime

figdir = "/home/ctroupin/DataOceano/MyOcean/figures/mooring/"
mooringdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/history/mooring/2plot/"
mooringlist = sorted(glob.glob(mooringdir+'*.nc'))

def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

# Prepare the new xticks (every month)
newxticks = []
monthname = []
for m in range(1, 13):
    newxticks.append((datetime.date(2015, m, 1) - datetime.date(2015, 1, 1)).days)
    monthname.append((datetime.date(2015, m, 1).strftime('%B')))

print monthname
print newxticks
print ' '

# Loop on the list of moorings
for moorings in mooringlist:
    print "working on "
    print moorings
    figname = os.path.basename(moorings).split('.')[0]+'_monthly'
    print figname

#   Load variables
    with netCDF4.Dataset(moorings) as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]
        ttime = nc.variables['TIME']
        timevec = ttime[:]
        time2plot = netCDF4.num2date(timevec, ttime.units)

#   Take only temperature with QC=1
        temp = nc.variables['TEMP'][:]
        tempQC = nc.variables['TEMP_QC'][:]
        temp = np.ma.masked_where(tempQC!=1, temp)

#    Load the dates with pandas
    ddates = pd.DatetimeIndex(time2plot)

#   Select days before and after 1st January 2015
    yearindex2015 = np.where(ddates.year == 2015)[0]
    yearindexpast = np.where(ddates.year < 2015)[0]

#   Convert the time to day of year to overlay the seasonal cycle
    day2plot = ddates.dayofyear
    month2plot = ddates.month

#   Initalize mask
    daymask = np.ma.make_mask_none(temp.shape)
    temp_masked = np.ma.array(temp)

#   Compute min, max and mean for each day of the year
    tempmin, tempmax, tempmean, tempstd = [], [], [], []
    maskedvalues = []

#   Loop on the months
    for mm in range(0, 12):

        print 'Working on month ' + monthname[mm]

#       Select values in the corresponding month, in the years before 2015
        monthindex = np.where(month2plot == mm+1)[0]
        monthindex = np.intersect1d(monthindex, yearindexpast)
        temp2calc = temp[monthindex]

#       Mask measurements outside 3 standards deviations
        values2mask = np.where(abs(temp2calc-temp2calc.mean())>=3*temp2calc.std())[0]
        daymask[monthindex[values2mask]] = True
        temp_masked = np.ma.masked_array(temp, daymask)

#       Compute min, max, mean and std for further use
        temp2calc = temp_masked[monthindex]
        tempmin.append(temp2calc.min())
        tempmax.append(temp2calc.max())
        tempmean.append(temp2calc.mean())
        tempstd.append(temp2calc.std())

#       Count how many bad values were removed
        maskedvalues.append(len(values2mask))

#   Put the mean and std in arrays
    tempmean = np.array(tempmean)
    tempstd = np.array(tempstd)
    tempmin = np.array(tempmin)
    tempmax = np.array(tempmax)

    tempmean2015 = []
#   Now work in 2015 (no need to loop on all the days)
    for mm in range(0, 8):
        monthindex = np.where(month2plot == mm+1)[0]
        monthindex = np.intersect1d(monthindex, yearindex2015)
        temp2calc = temp[monthindex]

        values2mask = np.where(abs(temp2calc-tempmean[mm])>=4*tempstd[mm])[0]
        daymask[monthindex[values2mask]] = True
        temp_masked = np.ma.masked_array(temp, daymask)

        temp2calc = temp_masked[monthindex]
        tempmean2015.append(temp2calc.mean())

    tempmean2015 = np.array(tempmean2015)


    tempmin = np.hstack((tempmin[-1], tempmin, tempmin[0]))
    tempmax = np.hstack((tempmax[-1], tempmax, tempmax[0]))
    tempmean = np.hstack((tempmean[-1], tempmean, tempmean[0]))
    #tempmean2015 = np.hstack((tempmean2015[-1], tempmean2015, tempmean2015[0]))
    time4plot = np.arange(-0.5, 13)

#   Make the plot
    fig = plt.figure(figsize=(14,7))
    ax = fig.add_subplot(111)

    plt.plot(time4plot, tempmin, 'x:', color='0.25', lw=2, label='2002-2014\nMinimum/Maximum')
    plt.plot(time4plot, tempmax, 'x:', color='0.25', lw=2)
    plt.plot(time4plot, tempmean, 'o-', color='0.25', lw=2, zorder=3, ms=5, label='2002-2014\nAverage')
    ax.fill_between(time4plot, tempmin, tempmax, facecolor='0.75', zorder=2, alpha=0.25)
    # #ax.fill_between(range(0, 366), tempmean-2*tempstd, tempmean+2*tempstd,facecolor='0.25', zorder=2, alpha=0.25)
    # #plt.plot(day2plot[yearindex], temp2[yearindex], 'ko', ms=1, markerfacecolor='0.5', lw=0.5)
    plt.plot(np.arange(0.5, len(tempmean2015)), tempmean2015, 'ro-',
             markerfacecolor='r', markeredgecolor='r', ms=7, lw=2, label='2015', zorder=4)


    plt.annotate('', xy=(5.5, tempmean[6]), xycoords='data', xytext=(5.5, tempmean2015[5]), textcoords='data',
                 arrowprops={'arrowstyle': '<->', 'color': '0.45', 'lw': 2})


    #plt.plot([5, 5], [tempmean[5], tempmean2015[5]], 'k--')
    text4plot = '$\Delta$T = ' + str(abs(np.round(tempmean[6]-tempmean2015[5], 1))) + '$^{\circ}$C'
    plt.text(4.7, 0.5*(tempmean[6]+tempmean2015[5]), text4plot, ha='right', color='0.45', fontsize=18)

    ax.set_xlim(0., 12.)
    ax.set_ylim(tempmin.min()-0.2, tempmax.max()+0.2)
    plt.xticks(np.arange(0, 12), monthname)
    fig.autofmt_xdate()
    plt.ylabel('Sea water \n temperature \n ($^{\circ}$C)', rotation=0, ha='right')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

#   Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.savefig(figdir + figname)
#    plt.show()
    plt.close()
