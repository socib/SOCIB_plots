#!/usr/bin/env python
#
# wind_plot_functions
# 
# --------------------------------------------------------------------------------------
def read_L2_wind_quikscat(windfile, domain):
    import numpy as np
    import netCDF4 as netcdf
    import time


    # Open NetCDF file
    with netcdf.Dataset(windfile) as nc:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        windspeed = nc.variables['retrieved_wind_speed'][:]
        winddirection = nc.variables['retrieved_wind_direction'][:]
        windtime = nc.variables['time'][:]

    # Change longitude
    lon[lon > 180] = lon[lon > 180] - 360.0

    # Reduce dimensions
    lon = lon.flatten()
    lat = lat.flatten()
    windspeed = windspeed.flatten()
    winddirection = winddirection.flatten()

    # ------------------------------------------------------------------------------------
    # Select sub-region
    goodlon = np.nonzero(np.logical_and(lon <= domain[1] + domain[4], lon >= domain[0] - domain[4]))
    goodlon = goodlon[0]

    if goodlon.size != 0:
        lat = lat[goodlon]
        lon = lon[goodlon]
        windspeed = windspeed[goodlon]
        winddirection = -winddirection[goodlon] + 90
        goodlat = np.nonzero(np.logical_and(lat <= domain[3] + domain[5], lat >= domain[2] - domain[5]))
        goodlat = goodlat[0]
        if goodlat.size != 0:
            lat = lat[goodlat]
            lon = lon[goodlat]
            windspeed = windspeed[goodlat]
            winddirection = winddirection[goodlat]

            uwind = windspeed * np.cos(np.deg2rad(winddirection))
            vwind = windspeed * np.sin(np.deg2rad(winddirection))
            uwind = np.ma.masked_where(uwind == uwind.fill_value, uwind)
            vwind = np.ma.masked_where(vwind == vwind.fill_value, vwind)
            uwind.data[uwind.data == uwind.data.min()] = 0
            vwind.data[vwind.data == vwind.data.min()] = 0


        else:

            print "No value in selected region"
            print " "
            lon, lat, uwind, vwind = [], [], [], []
    else:
        print "No value in selected region"
        print " "
        lon, lat, uwind, vwind = [], [], [], []

    return lon, lat, uwind, vwind, windtime


def plot_wind_sat(windfiles, domain, lon, lat, uwind, vwind, windtime, figdir, m):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as dt
    import time

    from matplotlib import colors


    figname = os.path.basename(windfiles)[:-3].replace('.', '_')

    cdict = {'red': ((0.0, 1.0, 1.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),
                       (0.5, .0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.20, 0.2),
                      (0.5, 0.0, 0.0),
                      (1.0, 0.0, 0.0))}

    cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    cmapvort = plt.cm.RdBu_r

    vmin = 0.
    vmax = 15.
    dvar = 2.5

    deltatime = (dt.datestr2num('1990,1,1') - dt.datestr2num('1970,1,1')) * 86400
    figtime = time.gmtime(windtime[0] + deltatime)


    # figtitle=str(figtime.tm_year)+'-'+str(figtime.tm_mon)+'-'+str(figtime.tm_mday)
    # figtitle = time.strftime("%d-%m-%Y %H:%M", time.gmtime(windtime[0, 0] + deltatime))

    bounds = np.arange(vmin, vmax + 0.01, dvar)
    bounds2 = np.arange(vmin + 0.1, vmax + 0.01, dvar)
    norm = colors.Normalize(vmin=vmin, vmax=vmax + 0.001)
  #  normvort = colors.Normalize(vmin=vmin, vmax=vmax + 0.001)


    meridians = np.arange(np.ceil(domain[0]), domain[1], domain[4])
    parallels = np.arange(domain[2], domain[3], domain[5])

    # Wind norm (should be same as windspeed)
    windnorm = np.sqrt(uwind * uwind + vwind * vwind)

    x, y = m(lon, lat)

    # Make the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m.ax = ax

    # contour=m.contourf(x,y,windspeed,levels2plot,cmap=cmap,norm=norm,extend='max')

    Q = m.quiver(x, y, uwind, vwind, windnorm,
                 units='width', scale=600, width=0.0015, norm=norm, cmap=cmap, zorder=3)

    ##    k = plt.quiverkey(Q, .9, 1.05, 10, r'$10\, ms^{-1}$',labelpos='W',
    ##                   fontproperties={'weight':'bold','size':'16'},color='k')

    m.drawparallels(parallels, linewidth=0.0, labels=[1, 0, 0, 0], fontsize=16)
    m.drawmeridians(meridians, linewidth=0.0, labels=[0, 0, 0, 1], fontsize=16)

    m.fillcontinents(ax=ax, color='0.5', zorder=2)
    m.drawcoastlines(color="k", ax=ax, linewidth=0.2, zorder=2)

    # Add the colorbar
    cbar = fig.colorbar(Q, cmap=cmap, orientation='vertical',
                        pad=0.025, aspect=15, shrink=1, norm=norm, extend='max')
    ax.annotate(r'($ms^{-1}$)', xy=(0.98, 1.03), xycoords='axes fraction')
    cbar.set_ticks(bounds)
    cbar.ax.set_yticklabels(bounds, fontsize=16)

    # plt.title(figtitle, fontproperties=prop, fontsize=24)
    ## Export figure and display it
    plt.savefig(figdir + figname, dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()


    # Make the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m.ax = ax

    # contour=m.contourf(x,y,windspeed,levels2plot,cmap=cmap,norm=norm,extend='max')

    Q = m.quiver(x, y, uwind, vwind, vort,
                 units='width', scale=600, width=0.0015, cmap=cmapvort, zorder=3)

    m.drawparallels(parallels, linewidth=0.0, labels=[1, 0, 0, 0], fontsize=16)
    m.drawmeridians(meridians, linewidth=0.0, labels=[0, 0, 0, 1], fontsize=16)

    m.fillcontinents(ax=ax, color='0.5', zorder=2)
    m.drawcoastlines(color="k", ax=ax, linewidth=0.2, zorder=2)

    # Add the colorbar
    cbar = fig.colorbar(Q, cmap=cmapvort, orientation='vertical',
                        pad=0.025, aspect=15, shrink=1, norm=norm, extend='both')
    ax.annotate(r'($ms^{-1}$)', xy=(0.98, 1.03), xycoords='axes fraction')
    cbar.set_ticks(bounds)
    cbar.ax.set_yticklabels(bounds, fontsize=16)

    # plt.title(figtitle, fontproperties=prop, fontsize=24)
    ## Export figure and display it
    plt.savefig(figdir + figname + '_curl', dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()


def plot_wind_sat_interp(windfiles, domain, lon, lat, uwind, vwind, windtime, figdir, m):

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as dt
    import time

    from matplotlib  import colors
    from mpl_toolkits.basemap import Basemap

    figname = os.path.basename(windfiles)[:-3].replace('.','_')

    cdict = {'red': ((0.0, 1.0, 1.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 1.0, 1.0),
                   (0.5, .0, 0.0),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.20, 0.2),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}

    cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
    cmapvort = plt.cm.RdBu_r
    vmin = 0.
    vmax = 15.
    dvar=2.5

    deltatime=(dt.datestr2num('1990,1,1')-dt.datestr2num('1970,1,1'))*86400

    bounds = np.arange(vmin,vmax+0.01,dvar)
    norm = colors.Normalize(vmin=vmin,vmax=vmax+0.001)
    normvort = colors.Normalize(vmin=-1.,vmax=1.0)
    meridians=np.arange(np.ceil(domain[0]),domain[1],domain[4])
    parallels=np.arange(domain[2],domain[3],domain[5])

    m = Basemap(projection='merc',llcrnrlon=domain[0],llcrnrlat=domain[2]+0.1,
                urcrnrlon=domain[1],urcrnrlat=domain[3],
                lat_ts=0.5*(domain[2]+domain[3]),
                resolution='l')

    # Wind norm (should be same as windspeed)
    windnorm=np.sqrt(uwind*uwind+vwind*vwind)

    x,y=m(lon,lat)

    # Make the plot
    fig=plt.figure()
    ax = fig.add_subplot(111)
    m.ax=ax

    # Wind curl
    f = 4*(np.pi/86400)*np.sin(np.deg2rad(0.5*(domain[2] * domain[3])))

    dx = np.gradient(x)
    dy = np.gradient(y)
    dux, duy = np.gradient(uwind)
    dvx, dvy = np.gradient(vwind)
    vort = (dvy / dx - dux / dy) / f

    Q=m.quiver(x, y, uwind, vwind, vort,
               units='width', scale=400, width=0.0015, norm=normvort, cmap=cmapvort)

    m.drawparallels(parallels, linewidth=0.0, labels=[1, 0, 0, 0], fontsize=16)
    m.drawmeridians(meridians, linewidth=0.0, labels=[0, 0, 0, 1], fontsize=16)

    m.fillcontinents(color="w", ax=ax, zorder=2)

    # Add the colorbar
    cbar=fig.colorbar(Q,cmap=cmap, orientation='vertical',
                      pad=0.025, aspect=15, shrink=1., norm=norm, extend='both')
    ax.annotate(r'($ms^{-1}$)', xy=(0.98, 1.03), xycoords='axes fraction')
    cbar.set_ticks(bounds)
    cbar.ax.set_yticklabels(bounds,fontsize=16)


    ## Export figure and display it
    plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)
    #plt.show()
    plt.close()