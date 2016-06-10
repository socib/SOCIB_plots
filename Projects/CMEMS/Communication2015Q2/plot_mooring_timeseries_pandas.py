__author__ = 'ctroupin'

import glob
import os
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime

plotmpl, plotpd, plotall = 0, 0, 1


figdir = "/home/ctroupin/DataOceano/MyOcean/figures/mooring/"
figdir2 = "/home/ctroupin/Projects/201501_InsTAC/MaterialQ2/"
mooringdir = "/home/ctroupin/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/history/mooring/2plot/"
mooringlist = sorted(glob.glob(mooringdir+'*.nc'))

# Prepare the new xticks (every month)
newxticks = []
monthname = []
for m in range(1, 13):
    newxticks.append((datetime.date(2015, m, 1) - datetime.date(2015, 1, 1)).days)
    monthname.append((datetime.date(2015, m, 1).strftime('%B')))

# Loop on the list of moorings
for moorings in mooringlist[0:1]:
    print "Working on file"
    print moorings
    print ' '
    figname = os.path.basename(moorings).split('.')[0]+'_monthly'
    print figname

#   Load variables
    with netCDF4.Dataset(moorings) as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]
        ttime = nc.variables['TIME']
        timevec = ttime[:]
        time2plot = netCDF4.num2date(timevec, ttime.units).squeeze()

#   Take only temperature with QC=1
        temperature_units = nc.variables['TEMP'].units
        temperature = nc.variables['TEMP'][:]
        temperatureQC = nc.variables['TEMP_QC'][:]


temperature[temperatureQC!=1] = np.nan
data = np.array((temperature[:, 0], temperatureQC[:, 0]))

# Create the data frame
hs = pd.DataFrame(data.T, index=time2plot, columns=('Temperature', 'Temperature QF'))

# Compute the monthly mean, min and max
hsmonth = hs['Temperature'].resample('M', how='mean')
hstempmax = hs['Temperature'].resample('M', how='max')
hstempmin = hs['Temperature'].resample('M', how='min')

# Compute the mean cycle (one value per month)
tempmean = []
tempmean2 = []
tempmin = []
tempmax = []
for mm in range(0, 12):
    tempmean.append(hsmonth[(hsmonth.index.month == mm+1)].mean())
    tempmean2.append(hsmonth[(hsmonth.index.month == mm+1) & (hsmonth.index.year < 2015)].mean())
    tempmin.append(hsmonth[(hsmonth.index.month == mm+1) & (hsmonth.index.year < 2015)].min())
    tempmax.append(hsmonth[(hsmonth.index.month == mm+1) & (hsmonth.index.year < 2015)].max())

tempmin = np.array(tempmin)
tempmax = np.array(tempmax)
tempmean= np.array(tempmean)

tempmin = np.hstack((tempmin[-1], tempmin, tempmin[0]))
tempmax = np.hstack((tempmax[-1], tempmax, tempmax[0]))
tempmean = np.hstack((tempmean[-1], tempmean, tempmean[0]))
#tempmean2015 = np.hstack((tempmean2015[-1], tempmean2015, tempmean2015[0]))
time4plot = np.arange(-0.5, 13)

print time4plot.shape
print tempmin.shape

if plotall:

    fig = plt.figure(figsize=(14,7))
    ax = fig.add_subplot(111)

    plt.plot(time4plot, tempmin, 'x:', color='0.25', lw=2, label='2002-2014\nMinimum/Maximum')
    plt.plot(time4plot, tempmax, 'x:', color='0.25', lw=2)
    plt.plot(time4plot, tempmean, 'o-', color='0.25', lw=2, zorder=3, ms=5, label='2002-2014\nAverage')
    ax.fill_between(time4plot, tempmin, tempmax, facecolor='0.75', zorder=2, alpha=0.25)
    #ax.fill_between(range(0, 366), tempmean-2*tempstd, tempmean+2*tempstd,facecolor='0.25', zorder=2, alpha=0.25)
    # #plt.plot(day2plot[yearindex], temp2[yearindex], 'ko', ms=1, markerfacecolor='0.5', lw=0.5)
    plt.plot(np.arange(0.5, 8), hsmonth['2015'], 'ro-',
             markerfacecolor='r', markeredgecolor='r', ms=7, lw=2, label='2015', zorder=4)

    plt.annotate('', xy=(5.5, tempmean[6]), xycoords='data', xytext=(5.5, hsmonth['2015-06']), textcoords='data',
                 arrowprops={'arrowstyle': '<->', 'color': '0.45', 'lw': 2})
    plt.plot([5, 5], [hsmonth['2015-05'], hsmonth['2015-05']], 'k--')
    text4plot = '$\Delta$T = ' + str(abs(np.round(hsmonth['2015-06'].values-hsmonth['2015-05'].values, 1)[0])) + '$^{\circ}$C'

    plt.text(4.7, 0.5*(tempmean[6]+hsmonth['2015-06']), text4plot, ha='right', color='0.45', fontsize=22)

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
    plt.savefig(figdir2 + figname + '.tif')
#    plt.grid(zorder=1)
    plt.show()
    plt.close()

