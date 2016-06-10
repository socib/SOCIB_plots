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
mooringlist = sorted(glob.glob(mooringdir+'*61198.nc'))

# Prepare the new xticks (every month)
newxticks = []
monthname = []
for m in range(1, 13):
    newxticks.append((datetime.date(2015, m, 1) - datetime.date(2015, 1, 1)).days)
    monthname.append((datetime.date(2015, m, 1).strftime('%B')))


yearmin, yearmax = 1998, 2015

kk = 0

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

    tempmin, tempmax, tempmean = [], [], []
    timemonth = []
    # Loop on years
    for yy in range(yearmin, yearmax+1):
#       Select days before and after 1st January 2015
        yearindex = np.where(ddates.year == yy)[0]

#       Loop on the months
        for mm in range(0, 12):

            print 'Working on month ' + monthname[mm]

    #       Select values in the corresponding month, in the years before 2015
            monthindex = np.where(ddates.month == mm+1)[0]
            monthindex = np.intersect1d(monthindex, yearindex)

            if len(monthindex) == 0:
                print 'Nothing'
                tempmin.append(np.nan)
                tempmax.append(np.nan)
                tempmean.append(np.nan)

            else:
                temp2calc = temp[monthindex]

        #       Compute min, max, mean and std for further use
                tempmin.append(temp2calc.min())
                tempmax.append(temp2calc.max())
                tempmean.append(temp2calc.mean())

            timemonth.append(mm+1)

            kk+=1

#   Put the mean and std in arrays
    tempmean = np.array(tempmean)
    tempmax = np.array(tempmax)
    tempmin = np.array(tempmin)
    timemonth = np.array(timemonth)

fig = plt.figure(figsize=(15, 7))
plt.plot(timemonth, tempmean, 'ko-')
# plt.plot(tempmax, 'r')
# plt.plot(tempmin, 'b')
plt.show()
plt.close()

print tempmean