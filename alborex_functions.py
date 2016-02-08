import numpy as np
import netCDF4 as netcdf
import cf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
# from okean import gshhs

import six
from six.moves import xrange

import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png


def prepare_map(coordinates, res):
    m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)
    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    return fig, m, ax


def prepare_map0(coordinates, res):
    m = Basemap(llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3], resolution=res)
    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    return fig, m, ax


def add_logo_on_map(imagepath, ax, position, zoom, zorder):
    logo2plot = read_png(imagepath)
    imagebox = OffsetImage(logo2plot, zoom=zoom)
    # coordinates to position this image

    ab = AnnotationBbox(imagebox, position,
                        xybox=(0., 0.),
                        xycoords='data',
                        pad=0.0,
                        boxcoords="offset points")
    ab.zorder = zorder

    ax.add_artist(ab)


def create_rect_patch(coordinates, m, alpha):
    xr1, yr1 = m(coordinates[0], coordinates[2])
    xr2, yr2 = m(coordinates[0], coordinates[3])
    xr3, yr3 = m(coordinates[1], coordinates[3])
    xr4, yr4 = m(coordinates[1], coordinates[2])
    verts = [(xr1, yr1), (xr2, yr2), (xr3, yr3), (xr4, yr4), (xr1, yr1), ]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY, ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='k', lw=0.1, alpha=alpha)
    return patch


def extract_coastline(coordinates, outputfile, coast_res):
    x, y = gshhs.get_coastline(xlim=[coordinates[0], coordinates[1]], ylim=[coordinates[2], coordinates[3]],
                               res=coast_res)
    # Save the coastline
    np.savetxt(outputfile, np.ma.vstack((x, y)).T)


def alborex_load_coast(coastdir, coastfile, valex):
    lon, lat = np.loadtxt(coastdir + coastfile, usecols=(0, 1), unpack=True)
    lon[lon == valex] = np.nan
    lat[lat == valex] = np.nan
    return lon, lat


def alborex_load_coast_gshhs(coastfile, coordinates, valex):
    lon, lat = np.loadtxt(coastfile, usecols=(0, 1), unpack=True)
    goodcoord = np.where((lon >= coordinates[0]) & (lon <= coordinates[1]) &
                         (lat >= coordinates[2]) & (lat <= coordinates[3]))[0]
    goodcoord2 = np.where(np.logical_or((lon == valex), (lat == valex)))[0]
    goodcoord = np.union1d(goodcoord, goodcoord2)
    lon, lat = lon[goodcoord], lat[goodcoord]
    lon[lon == valex] = np.nan
    lat[lat == valex] = np.nan
    return lon, lat


def alborex_load_bathy(bathydir, bathyfile, coordinates):
    with netcdf.Dataset(bathydir + bathyfile, 'r') as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        depth = depth[goodlat, :]
        depth = depth[:, goodlon]
    return lon, lat, depth


def alborex_load_altimetry(altimetryfile, coordinates):
    with netcdf.Dataset(altimetryfile, 'r') as nc:
        lon = nc.variables['lon'][:] - 360.
        lat = nc.variables['lat'][:]
        u = np.squeeze(nc.variables['u'][:])
        v = np.squeeze(nc.variables['v'][:])
        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        u = u[goodlat, :]
        u = u[:, goodlon]
        v = v[goodlat, :]
        v = v[:, goodlon]
    return lon, lat, u, v


def alborex_load_altimetry_time(altimetryfile, coordinates):
    with netcdf.Dataset(altimetryfile, 'r') as nc:
        time0 = nc.variables['time']
        time = netcdf.num2date(time0[:], time0.units)[0]
        lon = nc.variables['lon'][:] - 360.
        lat = nc.variables['lat'][:]
        u = np.squeeze(nc.variables['u'][:])
        v = np.squeeze(nc.variables['v'][:])
        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        u = u[goodlat, :]
        u = u[:, goodlon]
        v = v[goodlat, :]
        v = v[:, goodlon]
    return lon, lat, u, v, time


def alborex_load_sst(sstfile, coordinates):
    with netcdf.Dataset(sstfile, 'r') as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        timesst = nc.variables['time'][:]

        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        sst = np.squeeze(nc.variables['mcsst'][:, goodlat, goodlon])
        mask = nc.variables['lsmask'][:, goodlat, goodlon].squeeze()
        sst = np.ma.masked_where(mask == 1, sst)

        timesst *= 60.
        sat = nc.satellite
        sensor = nc.sensor_name
        sensorsat = sensor.upper() + ' on ' + sat.upper()
    return lon, lat, sst, timesst, sensorsat


def alborex_load_sst_meteofrance(sstfile, coordinates):
    with netcdf.Dataset(sstfile, 'r') as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]

        # subset
        goodlon = np.where(np.logical_and((lon >= coordinates[0]), (lon <= coordinates[1])))[0]
        goodlat = np.where(np.logical_and((lat >= coordinates[2]), (lat <= coordinates[3])))[0]
        lon = lon[goodlon]
        lat = lat[goodlat]
        sst = np.squeeze(nc.variables['sea_surface_temperature'][:, goodlat, goodlon])
        sst = sst - 273.15       # Kelvin to Celsius
        lon, lat = np.meshgrid(lon, lat)

        return lon, lat, sst



def alborex_load_sst_L2(sstdir, sstfile):
    if 'SST4' in sstfile:
        sstname = 'Geophysical_Data_sst4'
    else:
        sstname = 'Geophysical_Data_sst'

    with netcdf.Dataset(sstdir + sstfile, 'r') as nc:
        lon = nc.variables['Navigation_Data_longitude'][:]
        lat = nc.variables['Navigation_Data_latitude'][:]
        sst = nc.variables[sstname][:] * 0.005

    return lon, lat, sst


def alborex_load_sst_flag_L2(sstfile):
    if 'SST4' in sstfile:
        sstname = 'Geophysical_Data_sst4'
        sstflagname = 'Geophysical_Data_qual_sst4'
    else:
        sstname = 'Geophysical_Data_sst'
        sstflagname = 'Geophysical_Data_qual_sst'

    with netcdf.Dataset(sstfile, 'r') as nc:
        lon = nc.variables['Navigation_Data_longitude'][:]
        lat = nc.variables['Navigation_Data_latitude'][:]
        sst = nc.variables[sstname][:] * 0.005
        sstflag = nc.variables[sstflagname][:]
        sst = np.ma.masked_where(sstflag > 1, sst)
        sstyear = nc.Start_Year
        sstday = nc.Start_Day
    return lon, lat, sst, sstflag, sstyear, sstday


def alborex_load_ctd(ctdfile):
    with netcdf.Dataset(ctdfile, 'r') as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        depth = nc.variables['DEPTH'][:]
        temp = nc.variables['WTR_TEM_01'][:]
        chloro = nc.variables['CHLO'][:]
        chloro = np.ma.masked_where(np.isnan(chloro), chloro)
    return lon, lat, depth, temp, chloro


def alborex_load_ctd_TS(ctdfile):
    with netcdf.Dataset(ctdfile, 'r') as nc:
        depth = nc.variables['DEPTH'][:]
        temp = nc.variables['WTR_TEM_01'][:]
        salt = nc.variables['SALT_01'][:]
    return depth, temp, salt


def alborex_load_thermosal(thermosalfile):
    with netcdf.Dataset(thermosalfile, 'r') as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        time0 = nc.variables['time'][:]
        temp = nc.variables['WTR_TEM_REM'][:]
        salt = nc.variables['SALT'][:]
        fluor = nc.variables['FLUOR'][:]
        fluor = np.ma.masked_where(np.isnan(fluor), fluor)
    return lon, lat, time0, temp, salt, fluor


def alborex_load_ctd_time(ctdfile):
    with netcdf.Dataset(ctdfile, 'r') as nc:
        lon = nc.variables['LON'][:]
        lat = nc.variables['LAT'][:]
        depth = nc.variables['DEPTH'][:]
        temp = nc.variables['WTR_TEM_01'][:]
        chloro = nc.variables['CHLO'][:]
        chloro = np.ma.masked_where(np.isnan(chloro), chloro)
        time = nc.variables['time'][:]
    return lon, lat, depth, temp, chloro, time


def alborex_loadglider(gliderfile):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        temperature = nc.variables['temperature'][:]
    return lon, lat, depth, temperature


def alborex_loadglider_time(gliderfile):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        time = nc.variables['time'][:]
    return lon, lat, time


def alborex_loadglider_subsample(gliderfile, NN):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        lon = nc.variables['longitude'][::NN]
        lat = nc.variables['latitude'][::NN]
        depth = nc.variables['depth'][::NN]
        temperature = nc.variables['temperature'][::NN]
    return lon, lat, depth, temperature


def alborex_loadglider_coord(gliderfile):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
    return lon, lat, depth


def alborex_loadglider_TSP(gliderfile):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        salinity = nc.variables['salinity_corrected_thermal'][:]
        temperature = nc.variables['temperature'][:]
        pressure = nc.variables['pressure'][:]
        chloro = nc.variables['chlorophyll'][:]
    return temperature, salinity, pressure, chloro

def alborex_loadglider_oxy(gliderfile):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        oxygen = nc.variables['oxygen_concentration'][:]
    return oxygen


def alborex_loadglider_varname(gliderfile, varname, NN):
    with netcdf.Dataset(gliderfile, 'r') as nc:
        var = nc.variables[varname][::NN]
    return var


def change_wall_prop(ax, coordinates, depths, angles):
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_xaxis.gridlines.set_linestyles(':')
    ax.w_yaxis.gridlines.set_linestyles(':')
    ax.w_zaxis.gridlines.set_linestyles(':')
    ax.view_init(angles[0], angles[1])
    ax.set_xlim(coordinates[0], coordinates[1])
    ax.set_ylim(coordinates[2], coordinates[3])
    ax.set_zlim(depths[0], depths[1])
    ax.set_zlabel('Depth (m)')

    ax.set_zticks(np.arange(depths[0], depths[1] + 10, depths[2]))
    ax.set_zticklabels(range(int(-depths[0]), -int(depths[1]) - 10, -int(depths[2])))


def read_L2_wind(windfile, coordinates):
    # Open NetCDF file
    with netcdf.Dataset(windfile) as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        windspeed = nc.variables['wind_speed'][:]
        winddirection = nc.variables['wind_dir'][:]
        windtime = nc.variables['time'][:]
        bs_distance = nc.variables['bs_distance'][:]

    # Change longitude
    lon[lon > 180] = lon[lon > 180] - 360.0

    # Reduce dimensions
    lon = lon.flatten()
    lat = lat.flatten()
    windspeed = windspeed.flatten()
    winddirection = winddirection.flatten()
    bs_distance = bs_distance.flatten()

    # ------------------------------------------------------------------------------------
    # Select sub-region
    goodlon = np.nonzero(np.logical_and(lon <= coordinates[1], lon >= coordinates[0]))
    goodlon = goodlon[0]

    if goodlon.size != 0:
        lat = lat[goodlon]
        lon = lon[goodlon]
        windspeed = windspeed[goodlon]
        bs_distance = bs_distance[goodlon]
        winddirection = -winddirection[goodlon] + 90
        goodlat = np.nonzero(np.logical_and(lat <= coordinates[3], lat >= coordinates[2]))
        goodlat = goodlat[0]
        if goodlat.size != 0:
            lat = lat[goodlat]
            lon = lon[goodlat]
            windspeed = windspeed[goodlat]
            winddirection = winddirection[goodlat]
            bs_distance = bs_distance[goodlat]

            uwind = windspeed * np.cos(np.deg2rad(winddirection))
            vwind = windspeed * np.sin(np.deg2rad(winddirection))
            uwind = np.ma.masked_where(uwind == uwind.fill_value, uwind)
            vwind = np.ma.masked_where(vwind == vwind.fill_value, vwind)
            uwind.data[uwind.data == uwind.data.min()] = 0
            vwind.data[vwind.data == vwind.data.min()] = 0
        else:
            # print 'No value in selected region'
            # print ' '
            lon, lat, uwind, vwind = [], [], [], []
    else:
        print "No value in selected region"
        print " "
        lon, lat, uwind, vwind = [], [], [], []

    return lon, lat, uwind, vwind, windtime, bs_distance
