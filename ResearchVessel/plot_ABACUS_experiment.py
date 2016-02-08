#!/usr/bin/env python
# -------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors
import datetime, time, calendar
import matplotlib.text as text
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from alborex_functions import *
import pysocibclient

doplotsst, doplotvectors, doplotstream, doplotdrifters, doplotinset = 1, 0, 1, 1, 1


# Regions of interest and box for the experiment
dlon, dlat = .0, 1.0
coordinates = np.array((-6.75, 3.001, 34.75, 40.))
dlon2, dlat2 = .25, .25
coordinates2 = np.array((-1, -0.25, 36.65, 37.25))
coordinates3 = np.array((-6, -5.25, 35.8, 36.2))
res = 'h'

# Files and directories
# coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'
coastdir = '/home/ctroupin/Publis/201502_Alborex/data/coastline/'
coastfile = 'coastline_cartex.dat'
coastfile = 'coastline_cartex_f3.txt'

bathydir = '/home/ctroupin/IMEDEA/Cartex2014/data/bathymetry/'
bathyfile = 'topo_gebco_medsea.nc'
ctdfile = 'http://thredds.socib.es/thredds/dodsC/research_vessel/ctd/socib_rv-scb_sbe9002/L1/2014/dep0007_socib-rv_scb-sbe9002_L1_2014-05-25.nc'
gliderfile1 = 'http://thredds.priv.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L2/2014/dep0012_ideep00_ime-sldeep000_L2_2014-05-25_data_dt.nc'
gliderfile2 = 'http://thredds.priv.socib.es/thredds/dodsC/auv/glider/icoast00-ime_slcost000/L2/2014/dep0005_icoast00_ime-slcost000_L2_2014-05-25_data_dt.nc'

sstdir = '/home/ctroupin/DataOceano/Satellite/MODIS/data/L2/Alborex/SST/NetCDF/'
# altimetryfile = '/home/ctroupin/DataOceano/AVISO/netcdf_gridded/Alborex/dt_med_allsat_madt_uv_20140516_20141010.nc'

# sstdir = '/home/ctroupin/DataOceano/Satellite/SST/OceanColor/NetCDF/'
altimetryfile = '/home/ctroupin/DataOceano/AVISO/MedSea/Aviso_gridded/Alborex/dt_med_allsat_madt_uv_20140523_20141010.nc'

figdir = '/home/ctroupin/Publis/201502_Alborex/figures/MainFigure/'

sstfilelist = sorted(glob.glob(sstdir + 'A2014143020000*.nc'))
nfiles = len(sstfilelist)
valex = 999.

# Paths to icons
gliderlogo = '/home/ctroupin/Presentations/figures4presentations/icons/glider.png'
driferlogo = '/home/ctroupin/Presentations/figures4presentations/icons/surface_drifter.png'

# Time interval for the drifters (here approx. one week)
tt = datetime.datetime(2014, 5, 25, 0, 0, 0)
tt_end = datetime.datetime(2014, 6, 4, 0, 0, 0)
time_min = int(tt.strftime('%s'))
time_max = int(tt_end.strftime('%s'))
time_init = '2014-05-25T000000'
time_end = '2014-05-26T000000'

# Colormap
cmapsst = plt.cm.RdYlBu_r

# Compute min and max values 
sstmin, sstmax = 16., 21.
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .001, 1.0)

# Load coast
loncoast, latcoast = alborex_load_coast_gshhs(coastdir + coastfile, coordinates, valex)
#
# Load bathymetry
# lonbathy,latbathy,bathy = alborex_load_bathy(bathydir,bathyfile,coordinates)
#
# Load CTD and glider profiles
lonCTD, latCTD, depthCTD, tempCTD, sstCTD = alborex_load_ctd(ctdfile)
print "Loading CTD data"
print "..."

longlider1, latglider1, depthglider1 = alborex_loadglider_coord(gliderfile1)
print "Loading deep glider data"
print "..."
longlider2, latglider2, depthglider1 = alborex_loadglider_coord(gliderfile2)
print "Loading coastal glider data"
print "..."

# Generate lists of platforms
if doplotdrifters:
    socib_api = pysocibclient.ApiClient()
    drifterlist = socib_api.list_platforms(init_datetime=time_init, end_datetime=time_end,
                                           instrument_type="surface_drifter")

    print "Generating list of drifter files"
    print "..."

fig, m, ax = prepare_map(coordinates, res)
loncoast2, latcoast2 = m(loncoast, latcoast)
lonCTD, latCTD = m(lonCTD, latCTD)
longlider1, latglider1 = m(longlider1, latglider1)
longlider2, latglider2 = m(longlider2, latglider2)
plt.close()

for f in range(0, nfiles):

    fig, m, ax = prepare_map(coordinates, res)

    figname = "AlborexMission_V19"

    print figname
    # Load L2 data

    if doplotsst:
        lonSST, latSST, sst, sstflag, sstyear, sstday = alborex_load_sst_flag_L2(sstfilelist[f])
        print "Loading SST data"
        print "..."
        dd = datetime.datetime.strptime(str(sstyear) + ' ' + str(sstday), '%Y %j')
        sstdate = dd.strftime('%Y-%m-%d')

        x, y = m(lonSST, latSST)
        sstpcm = m.pcolormesh(x, y, sst, cmap=cmapsst, norm=normsst, edgecolor='none')
        cbar = fig.colorbar(sstpcm, cmap=cmapsst, norm=normsst, orientation='vertical', pad=-0.11, aspect=15,
                            shrink=0.6, extend='both')
        cbar.set_ticks(boundsst)
        cbar.set_label(r'$^{\circ}$C', fontsize=16, rotation=0)

    # Add currents from altimetry
    # ---------------------------

    lon_alti, lat_alti, u_alti, v_alti = alborex_load_altimetry(altimetryfile, coordinates)
    lon_alti, lat_alti = np.meshgrid(lon_alti, lat_alti)
    lon_alti, lat_alti = m(lon_alti, lat_alti)

    # Usual quiver
    if doplotvectors:
        Q1 = plt.quiver(lon_alti, lat_alti, u_alti, v_alti, width=0.001, scale=15, alpha=.5)
        qk = plt.quiverkey(Q1, 0.825, .05, 0.5, r'$0.5 m s^{-1}$', labelpos='E',
                           fontproperties={'weight': 'bold'})

    # Or streamfunction
    if doplotstream:
        speed = np.sqrt(u_alti.data * u_alti.data + v_alti.data * v_alti.data)
        speed = np.ma.masked_greater(speed, 1.5)
        plt.streamplot(lon_alti, lat_alti, u_alti, v_alti, color=speed,
                       arrowsize=3,
                       density=3, linewidth=.5, cmap=plt.cm.gray_r)

    # x and y ticks (no line plotted for better visibility)
    m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                    labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmapscale(-0.45, 35.1, -0.25, 35.1, 100, barstyle='simple', units='km', fontsize=12, zorder=3)

    # Labels
    xtext, ytext = m(coordinates[1] - .1, coordinates[2] + .25)
    #  plt.text(xtext, ytext, "AFRICA", va='center', ha='right',
    # fontsize=18, bbox=dict(fc='w', ec="w", lw=0.5,alpha=0.8, pad=5.0))

    xa1, ya1 = m(-5.7, 36.)
    xa2, ya2 = m(-6.1, 34.9)
    #    ax.annotate("Strait of\nGibraltar", xy=(xa1,ya1), xytext=(xa2,ya2),
    #            xycoords='data', textcoords='data',
    #            arrowprops=dict(arrowstyle="->",color='w',lw=10,
    #                            connectionstyle="arc3,rad=0.6",zorder=4),
    #            )
    #    ax.annotate("Strait of\nGibraltar", xy=(xa1,ya1), xytext=(xa2,ya2),
    #            xycoords='data', textcoords='data', fontsize=12,
    #            arrowprops=dict(arrowstyle="-",color='k',
    #                            connectionstyle="arc3,rad=0.",zorder=5),
    #            )

    # Add drifter trajectories on map
    # -------------------------------
    if doplotdrifters:

        londriftertotal, latdriftertotal, ndrifters = 0, 0, 0

        for drifter in drifterlist:
            drifter_opendap = drifter.product_list[-1].opendap
            if "2014-05-25" in drifter_opendap:
                with netcdf.Dataset(drifter_opendap) as nc:
                    timedrifter = nc.variables['time'][:]
                    goodtime = np.where((timedrifter > time_min) & (timedrifter < time_max))[0]
                    londrifter = nc.variables['LON'][goodtime]
                    latdrifter = nc.variables['LAT'][goodtime]
                    londrifter, latdrifter = m(londrifter, latdrifter)
                    m.plot(londrifter, latdrifter, 'ko', ms=0.1, alpha=.85, zorder=5)
                    londriftertotal += londrifter[0]
                    latdriftertotal += latdrifter[0]
                    ndrifters += 1
        add_logo_on_map(driferlogo, ax, [londriftertotal / ndrifters, latdriftertotal / ndrifters], 0.1, zorder=6)

    # Add rectangle for experiment
    patch = create_rect_patch(coordinates2, m, 0.1)
    ax.add_patch(patch)
    # and for Gibraltar
    patch = create_rect_patch(coordinates3, m, 0.2)
    ax.add_patch(patch)

    # CTD and gliders
    # m.plot(lonCTD,latCTD,'ko',ms=1, zorder=6)

    # Coastline and continent
    m.plot(loncoast2, latcoast2, 'k-', lw=0.1, zorder=4)
    m.fillcontinents(ax=ax, color='0.9', zorder=2)

    plt.title(r"25$-$31 May, 2014", fontsize=20)

    # Plot inset with CTD and gliders
    if doplotinset:
        axins = zoomed_inset_axes(ax, 4.5, loc=2)
        n1, n2 = m(coordinates2[0], coordinates2[2])
        n3, n4 = m(coordinates2[1], coordinates2[3])
        axins.set_xlim(n1, n3)
        axins.set_ylim(n2, n4)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

        # Add gliders tracks and CTD casts
        NN = 34

        axins.plot(lonCTD[:NN], latCTD[:NN], 'ks-', lw=0.5, color=".75", ms=5, zorder=2,
                   label=r"CTD 1st leg ($\times$ 34)", alpha=.7)
        axins.plot(lonCTD[NN:], latCTD[NN:], 'kx-', lw=0.2, ms=2.5, zorder=2, label=r"CTD 2nd leg ($\times$ 28)",
                   alpha=.7)

        axins.plot(longlider1, latglider1, '-', lw=3, color='k', label="Deep glider", zorder=3)
        axins.plot(longlider1, latglider1, '--', lw=3, color='.75', zorder=3)

        axins.plot(longlider2, latglider2, '-', lw=1.5, color='k', label="Coastal glider", zorder=4)
        axins.plot(longlider2, latglider2, '--', lw=1.5, color='.75', zorder=4)

        if doplotdrifters:
            axins.plot(londriftertotal / ndrifters, latdriftertotal / ndrifters, 'ko', ms=0.1,
                       label=r"Drifters ($\times$ 25)")
        if doplotsst:
            axins.pcolormesh(x, y, sst, cmap=cmapsst, norm=normsst, edgecolor='none')

        axins.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=12)

        add_logo_on_map(gliderlogo, axins, [longlider1[0], latglider1[0]], 0.1, zorder=6)
        add_logo_on_map(gliderlogo, axins, [longlider2[0], latglider2[0]], 0.1, zorder=6)

        mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")

    plt.savefig(figdir + figname, dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)

    # plt.show()
    plt.close()
