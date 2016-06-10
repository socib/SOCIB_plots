__author__ = 'ctroupin'

import glob
import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import datetime

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
    figname = os.path.basename(moorings).split('.')[0]
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

#   Select days before and after 1st January 2015
    day2015 = (datetime.datetime(2015,1,1, 0, 0) - datetime.datetime(1950,1,1, 0, 0)).total_seconds()/86400.
    yearindex2015 = np.where(timevec >= day2015)[0]
    yearindexpast = np.where(timevec < day2015)[0]

#   Convert the time to day of year to overlay the seasonal cycle
    day2plot = time2decimalday(time2plot)

#   Initalize mask
    daymask = np.ma.make_mask_none(temp.shape)
    temp_masked = np.ma.array(temp)

#   Compute min, max and mean for each day of the year
    tempmin, tempmax, tempmean, tempstd = [], [], [], []
    maskedvalues = []

    for dd in range(0, 366):

#       Select values in period of interest
        dayindex = np.where(np.logical_and((day2plot>=dd), (day2plot<dd+1)))[0]
        dayindex = np.intersect1d(dayindex, yearindexpast)
        temp2calc = temp[dayindex]

#       Mask measurements outside 3 standards deviations
        goodvalues = np.where(abs(temp2calc-temp2calc.mean())>=3*temp2calc.std())[0]
        daymask[dayindex[goodvalues]] = True
        temp_masked = np.ma.masked_array(temp, daymask)

#       Compute min, max, mean and std for further use
        temp2calc = temp_masked[dayindex]
        tempmin.append(temp2calc.min())
        tempmax.append(temp2calc.max())
        tempmean.append(temp2calc.mean())
        tempstd.append(temp2calc.std())

#       Count how many bad values were removed
        maskedvalues.append(len(goodvalues))

#   Put the mean and std in arrays
    tempmean = np.array(tempmean)
    tempstd = np.array(tempstd)

    tempmean2015 = []
#   Now work in 2015 (no need to loop on all the days)
    for dd in range(0, 250):

        dayindex = np.where(np.logical_and((day2plot>=dd), (day2plot<dd+1)))[0]
        dayindex = np.intersect1d(dayindex, yearindex2015)
        temp2calc = temp[dayindex]

        # Remove measurements outside 3 standards deviations
        goodvalues = np.where(abs(temp2calc-tempmean[dd])>=4*tempstd[dd])[0]
        daymask[dayindex[goodvalues]] = True
        temp_masked = np.ma.masked_array(temp, daymask)

        temp2calc = temp_masked[dayindex]
        tempmean2015.append(temp2calc.mean())

    tempmean2015 = np.array(tempmean2015)
    N = 5
#   Make the plot
    fig = plt.figure(figsize=(15, 6))
    ax = fig.add_subplot(111)

    plt.plot(range(0, 366), runningMeanFast(tempmin, N), ':', color='0.25', lw=2, label='2002-2014 Minimum/Maximum')
    plt.plot(range(0, 366), runningMeanFast(tempmax, N), ':', color='0.25', lw=2)
    plt.plot(range(0, 366), runningMeanFast(tempmean, N), '--', color='0.25', lw=2, zorder=3, label='2002-2014 Average')
    ax.fill_between(range(0, 366), runningMeanFast(tempmin, N), runningMeanFast(tempmax, N),
                    facecolor='0.75', zorder=2, alpha=0.25)
    #ax.fill_between(range(0, 366), tempmean-2*tempstd, tempmean+2*tempstd,facecolor='0.25', zorder=2, alpha=0.25)
    #plt.plot(day2plot[yearindex], temp2[yearindex], 'ko', ms=1, markerfacecolor='0.5', lw=0.5)
    plt.plot(range(0, 250), tempmean2015, 'ro-',
             markerfacecolor='r', markeredgecolor='r', ms=2, lw=2, label='2015', zorder=4)

    plt.legend(loc=2)
    ax.set_xlim(0., 366.-N)
    ax.set_ylim(10., 30.)
    plt.xticks(newxticks, monthname)
    fig.autofmt_xdate()
    plt.ylabel('Sea water \n temperature \n ($^{\circ}$C)', rotation=0, ha='right')
    plt.savefig(figdir + figname)
    plt.show()
    plt.close()
