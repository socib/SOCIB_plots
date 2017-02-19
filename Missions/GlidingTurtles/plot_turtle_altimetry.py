#!/usr/bin/env python
# -------------------------------------------------------------------------

from matplotlib import colors
import datetime
from alborex_functions import *

doplotsst, doplotstream, plotsalinity = 0, 1, 1

# Regions of interest and box for the experiment
dlon, dlat = 1.0, 1.0
coordinates = np.array((-1.00, 7.001, 36.9, 41.))
coordinates2 = np.array((-2.00, 8.001, 35.5, 42.))
res = 'f'
valex = -999.

figdate = "20150830"
figtitle = "30 August, 2015"

figname = ''.join(("Turtle_Altimetry_", figdate))
romsfile = ('http://thredds.priv.socib.es/thredds/dodsC/operational_models/'
            'oceanographical/hydrodynamics/wmop/latest.nc')
turtlefile = ('http://thredds.priv.socib.es/thredds/dodsC/animal/turtle/'
              'turtle_mel-alk_ttrk014/L1/2015/08/dep0001_turtle-mel_alk-ttrk014_L1_2015-08-12.nc')
altimetryfile = ('http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/'
                 'aviso/madt/L4/2015/08/nrt_med_allsat_madt_uv_{0}_{0}.nc.gz').format(figdate)
turtlelogo = "/home/ctroupin/Presentations/figures4presentations/logo/turtle.png"
figdir = '/home/ctroupin/Pictures/SOCIB/'

vmin, vmax, dvar = 36.5, 38.5, 0.5

variable = 'sea_surface_salinity'
cmap = plt.cm.RdYlBu_r

norm = colors.Normalize(vmin=vmin, vmax=vmax)
normvel = colors.Normalize(vmin=0.0, vmax=0.5)
bounds = np.arange(vmin, vmax + .0001, dvar)

# Compute min and max values 
sstmin, sstmax = 16., 21.
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .001, 1.0)

fig, m, ax = prepare_map(coordinates, res)
plt.close()

# Load turtle data
with netcdf.Dataset(turtlefile) as nc:
    lonturtle = nc.variables['LON'][:]
    latturtle = nc.variables['LAT'][:]
    lonturtleQC = nc.variables['QC_LON'][:]
    latturtleQC = nc.variables['QC_LAT'][:]

# Load model salinity
if plotsalinity:
    f = cf.read(romsfile)
    field2plot = f.select(variable)[0]
    lon = field2plot[0].coord('lon').array
    lat = field2plot[0].coord('lat').array
    field2plot = field2plot.array[0].squeeze()

    lon, lat = np.meshgrid(lon, lat)
    lon_ts, lat_ts = m(lon, lat)


fig, m, ax = prepare_map(coordinates, res)

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

lonturtle, latturtle = m(lonturtle, latturtle)

# Mask
lonturtle = np.ma.masked_where(lonturtleQC != 1, lonturtle)
latturtle = np.ma.masked_where(latturtleQC != 1, latturtle)

plt.plot(lonturtle, latturtle, 'ko-', zorder=5, markersize=1, linewidth=0.5)
plt.plot(lonturtle[-1], latturtle[-1], 'bo-', zorder=5, markersize=3, linewidth=0.5)

add_logo_on_map(turtlelogo, ax, [lonturtle[-1] , latturtle[-1]], 0.95, zorder=6)


# Add currents from altimetry
# ---------------------------

lon_alti, lat_alti, u_alti, v_alti = alborex_load_altimetry(altimetryfile, coordinates2)
lon_alti, lat_alti = np.meshgrid(lon_alti, lat_alti)
lon_alti, lat_alti = m(lon_alti, lat_alti)


# Or streamfunction
if doplotstream:
    speed = np.sqrt(u_alti.data * u_alti.data + v_alti.data * v_alti.data)
    speed = np.ma.masked_greater(speed, 1.5)
    plt.streamplot(lon_alti, lat_alti, u_alti, v_alti,
                   arrowsize=2,
                   density=7, linewidth=.25, color="0.5", zorder=4)

# x and y ticks (no line plotted for better visibility)
m.drawparallels(np.arange(round(coordinates[2]), coordinates[3], dlat), linewidth=0.,
               labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
m.drawmeridians(np.arange(round(coordinates[0]), coordinates[1], dlon), linewidth=0.,
                labels=[0, 0, 0, 1], fontname='Times New Roman', fontsize=16, zorder=1)

# m.drawmapscale(-0.45, 35.1, -0.25, 35.1, 100, barstyle='simple', units='km', fontsize=12, zorder=3)

pcm = m.pcolormesh(lon_ts, lat_ts, field2plot, norm=norm, cmap=cmap, zorder=3)

cbar = plt.colorbar(pcm, norm=norm, cmap=cmap, orientation='vertical', pad=0.05, aspect=15,
                        shrink=0.75, extend='both')
cbar.set_ticks(bounds)

# Coastline and continent
m.fillcontinents(ax=ax, color='0.9', zorder=2)
m.drawcoastlines(ax=ax, linewidth=0.1, color='k', zorder=4)
plt.title(figtitle, fontsize=20)


plt.savefig(figdir + figname, dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)

# plt.show()
plt.close()
