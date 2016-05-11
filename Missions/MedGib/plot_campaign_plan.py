import numpy as np
import glob
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib  import colors
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import splprep, splev

import medgib_functions

datadir = '/home/ctroupin/Publis/20150528_MedGib/data/'
datafiles = sorted(glob.glob(datadir+'*csv'))
coastdir = '/home/ctroupin/DataOceano/Coastlines/'
coastfile = 'medsea_coast2b.dat'
figdir = '/home/ctroupin/public_html/MedGib/'
figname = 'CampaignPlan'
valex = -999.

plotlarge, plotmedium, plotsmall = 0, 1, 0


coordinates1 = np.array((-6.25, -5., 35.75, 36.25))
dlon1,dlat1 = .25, .25

res1 = 'h'
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'
coastfile = 'coastline_cartex.dat'

sstmin,sstmax=16.,25.
cmapsst=plt.cm.RdYlBu_r
normsst= colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1,1.0)

# Prepare map and projection
fig, m, ax = medgib_functions.prepare_map(coordinates1, res1)

# Load coast
lonc, latc = medgib_functions.medgib_load_coast(coastdir, coastfile, valex)
lonc,latc=m(lonc,latc)
m.plot(lonc,latc,'k-',lw=0.5,zorder=4)

# Add SST measured by buoys
for f in datafiles:
    file2load=f
    code,year,month,day,hour,minute,second,lat,lon,sst,u,v,lat_qc,lon_qc,sst_qc,V_QC = np.loadtxt(file2load, delimiter=',', unpack=True)

    # Keep only good positions
    lon[lon_qc>1]=np.nan
    lon,lat=m(lon,lat)
    sstscat=plt.scatter(lon[0],lat[0],s=55,c=sst[0],edgecolor='none',cmap=cmapsst,norm=normsst)

m.drawparallels(np.arange(np.round_(coordinates1[2]),coordinates1[3],dlat1), linewidth=0.5,
                        labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16,zorder=1)
m.drawmeridians(np.arange(coordinates1[0],coordinates1[1],dlon1), linewidth=0.5,
                        labels=[0, 0, 1, 0], fontname='Times New Roman',fontsize=16,zorder=1)
m.drawmapscale(-4.,38.5,-4.,38,200, barstyle='simple', units='km', fontsize=12,zorder=4)


m.fillcontinents(ax=ax,color='w',zorder=3)
# plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close()
