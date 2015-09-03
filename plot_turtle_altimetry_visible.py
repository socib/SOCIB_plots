#!/usr/bin/env python
# -------------------------------------------------------------------------


from matplotlib import colors
import datetime
from alborex_functions import *
from osgeo import gdal

doplotsst, doplotstream, plotsalinity = 0, 1, 1

# Regions of interest and box for the experiment
dlon, dlat = 1.0, 1.0
coordinates = np.array((0.5, 6.50001, 36.9, 41.))
coordinates2 = np.array((-2.00, 8.001, 35.5, 42.))
valex = -999.

figdate = "20150902"
figtitle = "2 September, 2015"

figname = "Turtle_Altimetry_"+ figdate
romsfile = "http://thredds.priv.socib.es/thredds/dodsC/operational_models/oceanographical/hydrodynamics/wmop/latest.nc"
turtlefile = "http://thredds.priv.socib.es/thredds/dodsC/animal/turtle/turtle_mel-alk_ttrk014/L1/2015/08/dep0001_turtle-mel_alk-ttrk014_L1_2015-08-12.nc"
altimetryfile = "http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/madt/L4/2015/09/nrt_med_allsat_madt_uv_" + figdate + "_" + figdate + ".nc.gz"
altimetryfile2 = "http://thredds.priv.socib.es/thredds/dodsC/satellite/altimetry/aviso/madt/L4/2015/09/nrt_med_allsat_madt_h_" + figdate + "_" + figdate + ".nc.gz"
turtlelogo = "/home/ctroupin/Presentations/figures4presentations/logo/turtle.png"
figdir = '/home/ctroupin/Pictures/SOCIB/'

visiblefile = '/data_local/Satellite/Visible/nasa-worldview-2015-09-01.tiff'


gtif = gdal.Open(visiblefile)
gtif.GetProjectionRef()

#Set the plot axis limits to the proper map coordinates.
arr = gtif.ReadAsArray()
trans = gtif.GetGeoTransform()
extent = (trans[0], trans[0] + gtif.RasterXSize*trans[1],
          trans[3] + gtif.RasterYSize*trans[5], trans[3])

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

    lon_ts, lat_ts = np.meshgrid(lon, lat)

# Load altimetry
with netcdf.Dataset(altimetryfile) as nc:
    lon_alti = nc.variables['lon'][:] - 360.
    lat_alti = nc.variables['lat'][:]
    u = nc.variables['u'][:].squeeze()
    v = nc.variables['v'][:].squeeze()

with netcdf.Dataset(altimetryfile2) as nc:
    lon_alti2 = nc.variables['lon'][:] - 360.
    lat_alti2 = nc.variables['lat'][:]
    adt = nc.variables['adt'][:].squeeze()

adt = adt - adt.mean()
normalti = colors.Normalize(vmin=-0.15, vmax=0.15)


fig = plt.figure()
ax = plt.subplot(111)

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels

ax.set_xlim(coordinates[0], coordinates[1])
ax.set_ylim(coordinates[2], coordinates[3])

# Mask
# lonturtle = np.ma.masked_where(lonturtleQC >  2, lonturtle)
# latturtle = np.ma.masked_where(latturtleQC >  2, latturtle)

plt.imshow(arr[:3,:,:].transpose((1, 2, 0)), extent=extent, zorder=3,alpha=1)

plt.plot(lonturtle[3:], latturtle[3:], 'go-', zorder=5, markersize=1, linewidth=2.5)
add_logo_on_map(turtlelogo, ax, [lonturtle[-2] , latturtle[-2]], 1., zorder=6)


# Add currents from altimetry
# ---------------------------



# Or streamfunction
if doplotstream:

    plt.streamplot(lon_alti, lat_alti, u, v, color=adt, norm=normalti, arrowsize=2,
               density=12, linewidth=.25, cmap=plt.cm.RdBu_r, zorder=5)


plt.title(figtitle, fontsize=20)


plt.savefig(figdir + figname + '_visible', dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)

# plt.show()
plt.close()
