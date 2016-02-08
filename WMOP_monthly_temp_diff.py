#!/usr/bin/python
__author__ = 'ctroupin'

import netCDF4
import numpy as np
from calendar import monthrange
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import Basemap
from matplotlib import colors


def extract_lon_lat(datafile):
    with netCDF4.Dataset(datafile) as nc:
        lon = nc.variables['lon_rho'][:]
        lat = nc.variables['lat_rho'][:]
    return lon, lat


def compute_monthly_mean(year, month, path, nlon, nlat):
    temperature_sum = np.zeros((nlat, nlon))
    yy, mm = str(year), str(month).zfill(2)
    ndays = 0
    for day in range(1, monthrange(year, month)[-1] + 1):
        dd = str(day).zfill(2)
        datafile = "%s/%s/%s/roms_wmop_%s%s%s.nc" % (path, yy, mm, yy, mm, dd)
        print datafile

        try:
            with netCDF4.Dataset(datafile, 'r') as nc:
                temperature = nc.variables['temp'][:]
                np.ma.masked_outside(temperature, -1, 40., copy=True)
                temperature_sum += np.ma.average(temperature, 0)
                ndays += 1
        except RuntimeError:
            print("This file doesn't exist")

    return temperature_sum / ndays


def write_monthly_mean(field, newfile, nlon, nlat):
    rootgrp = netCDF4.Dataset(newfile, "w", format="NETCDF4")
    lat = rootgrp.createDimension("lat", nlat)
    lon = rootgrp.createDimension("lon", nlon)
    temp = rootgrp.createVariable("temp", "f4", ("lat", "lon",))
    temp[:] = field
    rootgrp.close()
    return


def make_plot_anomalies(coordinates, lon, lat, file1, file2, figname):
    with netCDF4.Dataset(file1) as nc:
        field1 = nc.variables['temp'][:]
    with netCDF4.Dataset(file2) as nc:
        field2 = nc.variables['temp'][:]

    fig, m, ax = prepare_map(coordinates, 'h')
    cmap = plt.cm.bwr
    norm = colors.Normalize(vmin=-2.5, vmax=2.5)
    anomalies = field1 - field2
    # np.ma.masked_where(anomalies==0.0, anomalies,copy=True)
    lon, lat = np.meshgrid(lon, lat)
    lon, lat = m(lon, lat)

    pcm = m.pcolormesh(lon, lat, field1 - field2, cmap=cmap, norm=norm)
    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], 2), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], 2), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)
    plt.colorbar(pcm, extend='both')
    m.fillcontinents(ax=ax, color='0.5', zorder=2)
    m.drawcoastlines(ax=ax, linewidth=0.2)
    plt.savefig(figname)
    plt.close()
    return


def prepare_map(coordinates, res):
    m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)
    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    return fig, m, ax


def main():

    year1, year2 = 2015, 2014
    month = 5
    coordinates = [-4, 10, 35, 44]
    path = "http://thredds.priv.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop"
    dir = "/home/ctroupin/SOCIB/Facilities/Modelling/figures/Anomalies/"
    romsfile = "http://thredds.priv.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop/2015/12/roms_wmop_20151231.nc"

    lon_roms, lat_roms = extract_lon_lat(romsfile)
    nlon, nlat = len(lon_roms), len(lat_roms)

    file1 = os.path.join(dir, 'temperature_avg_' + str(year1) + str(month).zfill(2) + '.nc')
    file2 = os.path.join(dir, 'temperature_avg_' + str(year2) + str(month).zfill(2) + '.nc')

    if not (os.path.exists(file1)):
        temperature_avg1 = compute_monthly_mean(year1, month, path, nlon, nlat)
        write_monthly_mean(temperature_avg1, file1, nlon, nlat)
    if not (os.path.exists(file2)):
        temperature_avg2 = compute_monthly_mean(year2, month, path, nlon, nlat)
        write_monthly_mean(temperature_avg2, file2, nlon, nlat)

    figname = "anomalies_%s%s_%s%s.png" % (str(year1), str(month), str(year2), str(month))
    figname = os.path.join(dir, figname)
    make_plot_anomalies(coordinates, lon_roms, lat_roms, file1, file2, figname)


if __name__ == '__main__':
    main()
