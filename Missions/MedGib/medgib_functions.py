import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def prepare_map(coordinates, res):
    m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2], 
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)
    fig = plt.figure()
    ax = plt.subplot(111)
    m.ax = ax
    return fig, m, ax


def medgib_load_coast(coastdir, coastfile, valex):
    lon, lat = np.loadtxt(coastdir + coastfile, usecols=(0, 1), unpack=True)
    lon[lon == valex] = np.nan
    lat[lat == valex] = np.nan
    return lon, lat


def medgib_load_bathy(bathydir, bathyfile, coordinates):
    with netCDF4.Dataset(bathydir + bathyfile, 'r') as nc:
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


def medgib_load_sst(sstfile, coordinates):
    with netCDF4.Dataset(sstfile, 'r') as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        timesst = nc.variables['time'][:]
        # lsmask = nc.variables['lsmask'][:]

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


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS = 0.00000000005
    if (np.abs(np.cos(glat1)) < EPS) and not (np.abs(np.sin(faz)) < EPS):
        print("Only N-S courses are meaningful, starting at a pole!")

    a = 6378.13 / 1.852
    f = 1 / 298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if cf == 0:
        b = 0.
    else:
        b = 2. * np.arctan2(tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while np.abs(y - c) > EPS:
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
             d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2 * np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2 * np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180. / np.pi
    glat2 *= 180. / np.pi
    baz *= 180. / np.pi

    return (glon2, glat2, baz)


def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    # ~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X, Y = m(X, Y)
    plt.plot(X, Y, **kwargs)
