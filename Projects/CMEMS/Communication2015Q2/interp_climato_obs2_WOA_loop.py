__author__ = 'ctroupin'


"""
We interpolate the WOA2013 (0.25 degree resolutio) at the location of the observations
then we compute the anomalies.
"""
import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import ndimage
from mpl_toolkits.basemap import Basemap
import calendar

plottest1, plottest2, plotanom = 0, 0, 1

# Select year, month, depth and platform
year = 2015
monthlist = [6, 7]
platformlist  = ['profiler-glider']
depthlist = [10., 20., 30., 50., 100.,]

figdir = "/home/ctroupin/DataOceano/MyOcean/figures/anomaliesWOA/V3/"
obsdir = "/home/ctroupin/DataOceano/MyOcean/4interp/V3/"

if not os.path.exists(figdir):
    os.makedirs(figdir)


cmap = plt.cm.hot_r
cmapanom = plt.cm.RdBu_r
fmin, fmax = 10., 30.
fmin2, fmax2 = 17., 28.
norm = colors.Normalize(vmin=fmin, vmax=fmax)
norm2 = colors.Normalize(vmin=fmin2, vmax=fmax2)
# normanom = colors.Normalize(vmin=-2.5, vmax=2.5)
boundsanom = np.arange(-4., 4.00001, 1.0)
coordinates = np.array((-6.75, 37., 29, 46.))
dlon, dlat = 5., 3.

normanom = colors.BoundaryNorm(boundsanom, cmapanom.N)

parallels2plot = range(30, 46, 5)
meridians2plot = range(-5, 36, 5)

mypad = 0.03
myaspect = 20

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], urcrnrlon=coordinates[1],
            urcrnrlat=coordinates[3], lat_ts=0.5*(coordinates[2]+coordinates[3]), resolution='f')


# Loop on platforms
for platform in platformlist:
    print 'Working on platform ' + platform
    # Loop on month
    for month in monthlist:
        print 'Working on month ' + str(month)

        period = str(year)+str(month).zfill(2)
        climatofile = "/data_local/DataOceano/WOA/med/woa13_decav_t" + str(month).zfill(2) + "_04v2.nc"

        # Loop on depth
        for depthinterp in depthlist:
            print 'Working on depth ' + str(depthinterp)


            outputfile = obsdir + 'anomalies_' + str(int(depthinterp)) + '_' + platform + '_' + period + '.dat'
            timefile = obsdir + 'anomaliestime_' + str(int(depthinterp)) + '_' + platform + '_' + period + '.dat'

            depthtext = 'Depth: ' + str(int(depthinterp)) + ' m'

            # Read the climato at the chosen depth
            with netCDF4.Dataset(climatofile) as nc:
                depthclim = nc.variables['depth'][:]
                lonclim = nc.variables['lon'][:]
                latclim = nc.variables['lat'][:]

                depthclimindex = np.where(depthclim == depthinterp)[0]
                tempclim2interp = nc.variables['t_an'][0, depthclimindex, :, :].squeeze()

            # Read the interpolated observations and their time (generated by interp_vessel_profiles)
            lonobs, latobs, tempobs = np.loadtxt(outputfile, unpack=True)
            yearobs, monthobs, daymonths = np.loadtxt(timefile, unpack=True)

            # print len(lonobs)
            # print len(yearobs)

            # Interpolate the climato on the positions of the observations
            # Grid necessary for interpolating

            llonclim, llatclim = np.meshgrid(lonclim, latclim)
            lonclim2 = llonclim.flatten(1)
            latclim2 = llatclim.flatten(1)

            # Interpolate climatology where we have coordinates using map_coordinates
            # (requires a small coordinates change)
            lonobs2 = (lonobs-lonclim[0])/(lonclim[1]-lonclim[0])
            latobs2 = (latobs-latclim[0])/(latclim[1]-latclim[0])

            # print tempclim2interp.shape
            # print lonobs2.shape
            # print latobs2.shape
            Climvalues = ndimage.map_coordinates(tempclim2interp.T, np.array((lonobs2, latobs2)), order=1, mode='nearest')


            # Make plots

            lonclim, latclim = m(llonclim, llatclim)
            lonobs, latobs = m(lonobs, latobs)

            if plottest1:

                figname = 'obs_climato_' + str(int(depthinterp)) +  platform + '_' + period
                fig = plt.figure()
                ax = fig.add_subplot(111)

                m.fillcontinents(color="0.4", lake_color="0.4", zorder=4)
                m.drawcoastlines(ax=ax, linewidth=0.2, color='0.35', zorder=5)

                m.drawparallels(parallels2plot, linewidth=0.,
                                labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
                m.drawmeridians(meridians2plot, linewidth=0.,
                                labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)
                scat = plt.scatter(lonobs, latobs, s=7, c=tempobs, cmap=cmap, norm=norm2, edgecolor='None', zorder=3)
                plt.pcolormesh(lonclim, latclim, tempclim2interp, cmap=cmap, norm=norm2, zorder=2)

                cbar = plt.colorbar(scat, orientation='horizontal', pad=mypad, aspect=20, extend='both', ax=ax)
                cbar.set_label('$^{\circ}$C')
                plt.text(0.02, 0.05, depthtext, color='w', fontsize=24, ha='left', va='center', transform=ax.transAxes, zorder=6)
                plt.text(0.02, 0.15, calendar.month_name[month], color='w', fontsize=24, ha='left', va='center', transform=ax.transAxes, zorder=6)
                plt.savefig(figdir + figname, dpi=300)

            #    plt.show()
                plt.close()

            if plottest2:

                figname = 'obs_interp_climato_' + str(int(depthinterp)) + platform + '_' + period

                fig = plt.figure()
                ax = fig.add_subplot(111)
                m.fillcontinents(color="0.4", lake_color="0.4", zorder=4)
                m.drawcoastlines(ax=ax, linewidth=0.2, color='0.35', zorder=5)

                m.drawparallels(parallels2plot, linewidth=0.,
                                labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
                m.drawmeridians(meridians2plot, linewidth=0.,
                                labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)

                plt.pcolormesh(lonclim, latclim, tempclim2interp, cmap=cmap, norm=norm2, zorder=2)
                scat = plt.scatter(lonobs, latobs, s=7, c=Climvalues, cmap=cmap, norm=norm2, edgecolor='None', zorder=3)
                plt.plot(lonobs, latobs, 'ko', ms=0.2, zorder=4)
                cbar = plt.colorbar(scat, orientation='horizontal', pad=mypad, aspect=myaspect, extend='both', ax=ax)
                cbar.set_label('$^{\circ}$C')
                plt.text(0.02, 0.05, depthtext, color='w', fontsize=24, ha='left', va='center', transform=ax.transAxes, zorder=6)
                plt.text(0.02, 0.15, calendar.month_name[month], color='w', fontsize=24, ha='left', va='center', transform=ax.transAxes, zorder=6)
                plt.savefig(figdir + figname, dpi=300)
            #    plt.show()
                plt.close()

            # Compute the anomalies
            tempanomalies = tempobs - Climvalues


            # Plot anomalies
            if plotanom:

                figname = 'anomalies_' + str(int(depthinterp)) + platform + '_' + period

                fig = plt.figure()
                ax = fig.add_subplot(111)
                m.fillcontinents(color="0.4", lake_color="0.4", zorder=4)
                m.drawcoastlines(ax=ax, linewidth=0.2, color='0.35', zorder=5)

                m.drawparallels(parallels2plot, linewidth=0.,
                                labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
                m.drawmeridians(meridians2plot, linewidth=0.,
                                labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)


                scat = plt.scatter(lonobs, latobs, s=7, c=tempanomalies, cmap=cmapanom, norm=normanom, edgecolor='None', zorder=3)
                cbar = plt.colorbar(scat, orientation='horizontal', pad=mypad, boundaries=boundsanom, aspect=myaspect, extend='both', ax=ax)
                cbar.set_label('Temperature anomaly ($^{\circ}$C) with respect to climatology')
                cbar.set_ticks(boundsanom)
                plt.text(0.02, 0.05, depthtext, color='w', fontsize=24, ha='left', va='center', transform=ax.transAxes, zorder=6)
                plt.text(0.02, 0.15, calendar.month_name[month] + ' ' + str(year), color='w', fontsize=24, ha='left', va='center', transform=ax.transAxes, zorder=6)
                plt.savefig(figdir + figname + '.png', dpi=300)
            #    plt.show()
                plt.close()


print 'Figures available in ' + figdir