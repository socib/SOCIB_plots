__author__ = 'ctroupin'

import glob
import netCDF4
import numpy as np
import sys

sys.path.append("/home/ctroupin/Software/Python/seawater-3.3.2")
import seawater

year = 2014
monthlist = [6, 7]
platformlist  = ['profiler-glider', 'vessel']
depthlist = [10., 20., 30., 50., 100.,]

outputdir = "/home/ctroupin/DataOceano/MyOcean/4interp/V3/"

# Loop on platforms
for platform in platformlist:
    print 'Working on platform ' + platform
    # Loop on month
    for month in monthlist:
        print 'Working on month ' + str(month)

        period = str(year)+str(month).zfill(2)
        datadir = "/data_local/DataOceano/MyOcean/INSITU_MED_NRT_OBSERVATIONS_013_035/monthly/" + platform + "/" + period + "/"
        filelist = sorted(glob.glob(datadir + '*.nc'))


        # Loop on depth
        for depthinterp in depthlist:

            print 'Working on depth ' + str(depthinterp)

            outputfile = outputdir + 'anomalies_' + str(int(depthinterp)) + '_' + platform + '_' + period  + '.dat'
            timefile = outputdir + 'anomaliestime_' + str(int(depthinterp)) + '_' + platform+ '_' + period + '.dat'

            loninterp, latinterp, tempinterp, timeinterp = [], [], [], []

            for datafile in filelist:
                print datafile

                with netCDF4.Dataset(datafile) as nc:
                    lon = nc.variables['LONGITUDE'][:]
                    lat = nc.variables['LATITUDE'][:]
                    time = nc.variables['TIME']
                    time2plot = netCDF4.num2date(time[:], time.units)
                    nprofiles = len(lon)
                    # print 'Number of profiles = ' + str(nprofiles)
                    try:
                        depth = nc.variables['DEPH'][:]
                    except KeyError:
                        # print "No variable depth"
                        try:
                            pressure = nc.variables['PRES'][:]
                            #print pressure.shape
                            #print "convert pressure to depth"
                            depth = np.zeros_like(pressure)
                            for ii in range(0, pressure.shape[0]):
                                depth[ii, :] = seawater.eos80.dpth(pressure[ii, :], lat[ii])
                        except KeyError:
                            # print "No variable pressure"
                            print ' '
                    ndepth = depth.shape[1]
                    if ndepth > 5:
                        print 'Number of depth levels = '+str(ndepth)
                        try:
                            temp = nc.variables['TEMP'][:]

                            for ii in range(0, nprofiles):
                            # Interpolate temperature at given level
                                print 'Depth min. = ' + str(depth[ii, :].min())
                                print 'Depth max. = ' + str(depth[ii, :].max())

                                if (depth[ii, :].min() <= depthinterp) & (depth[ii, :].max() >= depthinterp):
                                    # print "Ok to interpolate the profiles"
                                    gooddepth = np.where(depth[ii, :] > 0.)[0]
                                    tempinterp0 = np.interp(depthinterp, depth[ii, gooddepth], temp[ii, gooddepth])
                                    # print tempinterp0

                                    if ( (tempinterp0 < 35.) & (lon[ii] < 180.)):
                                        loninterp.append(lon[ii])
                                        latinterp.append(lat[ii])
                                        tempinterp.append(tempinterp0)
                                        timeinterp.append(time2plot[ii].strftime('%Y %m %d'))
                                    # fig = plt.figure()
                                    # ax = fig.add_subplot(111)
                                    # plt.plot(temp[ii, :], -depth[ii,:], 'k-', lw=1)
                                    # plt.plot(temp[ii,:].mean(), -depthinterp, 'ro', ms=15)
                                    # plt.plot(tempinterp0, -depthinterp, 'bo', ms=10)
                                    # plt.show()
                                    # plt.close()
                        except KeyError:
                            print "No variable depth"
                                #   Loop on the profiles



            #print timeinterp
            np.savetxt(outputfile, np.c_[loninterp, latinterp, tempinterp], fmt='%2.4f %2.5f %2.3f')
            np.savetxt(timefile, timeinterp, fmt='%s')


            # print 'Finished writing'
