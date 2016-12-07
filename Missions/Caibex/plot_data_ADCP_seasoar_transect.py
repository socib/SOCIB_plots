#--------------------------------------------------------------------------
#
# plot_data_ADCP_seasoar_transect.py
#
# plot the ADCP velocities along the Seasoar transect
# (no need to interpolate horizontally, maybe vertically)
#
# ctroupin, January 2015
#--------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
from scipy import interpolate
import matplotlib.dates as dt
import time, datetime
import scipy.io

os.system('clear')
basemap_resolution = 'f'
doplot = 1
frq = 75
vlimit = 0.8

#depthminlist = (10.,25.,50.,75.,100.)
#depthmaxlist = (10.,25.,50.,75.,100.)
depthminlist = (0.,)
depthmaxlist = (100.,)

# files and directories
seasoardir  = '/home/ctroupin/DataOceano/CAIBEX_campaign/data/Seasor/' 
datadir = '/home/ctroupin/DataOceano/CAIBEX_campaign/data/ADCP/' 
figdir = '/home/ctroupin/Publis/201406_CAIBEX/figures/2015/' 
figdir='/home/ctroupin/DataOceano/CAIBEX_campaign/data/ADCP/figures_python/SeaSoar_transect/'
figtype='.png'

if not(os.path.exists(figdir)):
    os.makedirs(figdir)

figbasename0 = 'adcp_seasoar'+str(frq)+'khz_V3' 

# --------------------------------------------------------------------------------------

# region of interest
lonmin = -12.25 
lonmax = -9.5 
latmin = 29.75
latmax = 31.75

# Where to put the lon/lat ticks
dlon,dlat=.5,.5
lon2plot = np.arange(-12.5,-9.,dlon)
lat2plot = np.arange(29.5,32.,dlon)

# Spacing between x/y labels
dlon = 0.5
dlat = 0.5

# prepare projection
m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                    urcrnrlon=lonmax,urcrnrlat=latmax,  \
                    lat_ts=0.5*(latmin+latmax),\
                    resolution=basemap_resolution)

# time for the Seasoar transect: 25-27 August 2009
timeref=dt.datestr2num('2009,1,1')
timestart1 = dt.datestr2num('2009,8,17')-timeref+1 
timeend1 = dt.datestr2num('2009,8,20')-timeref+1 


# --------------------------------------------------------------------------------------
# Load ADCP data
ADCPcoordfile=datadir+'CG'+str(frq)+'_xy.mat'
ADCPvaluefile=datadir+'CG'+str(frq)+'_uv.mat'
mat = scipy.io.loadmat(ADCPcoordfile)
xyt = mat['xyt']
zc = mat['zc']
mat = scipy.io.loadmat(ADCPvaluefile)
uv = mat['uv']

# Separate variables into vectors
u=uv[:,::2]
v=uv[:,1::2]    
lon=-360+xyt[0,:]
lat=xyt[1,:]
time=xyt[2,:]

# Mask nan values
u=np.ma.masked_array(u,np.isnan(u))
v=np.ma.masked_array(v,np.isnan(v))

# select good time
goodtime =  np.nonzero(np.logical_and(time>=timestart1,time<=timeend1))[0]
goodtime = range(0,932)    
lon1 = lon[goodtime] 
lat1 = lat[goodtime] 
u1 = u[:,goodtime] 
v1 = v[:,goodtime] 


### Mask largest values (unrealistic)
##goodval = np.nonzero(np.logical_and(abs(u1)<=vlimit,abs(v1)<=vlimit))[0]
##u1 = u1[goodval]
##v1 = v1[goodval]
##lon1 = lon1[goodval]
##lat1 = lat1[goodval]


# compute Coriolis frequency
f = 2*7.2921e-5*np.sin(np.deg2rad(30.5))

# Loop on a depth list

for depthmin, depthmax in zip(depthminlist,depthmaxlist):

    # find the good depth,
    if (depthmin == depthmax):
        titletext = 'SeaSoar tracks 1-6, depth: '+str(int(depthmin))+' m, '+str(frq)+' kHz'
        figbasename=figbasename0+'_SeaSoar'+str(int(depthmin))+'m'
        depth = depthmin
        gooddepth = (abs(zc-depth)).argmin()

        
    else:
        gooddepth = np.nonzero(np.logical_and(zc>=depthmin,zc<=depthmax))[0]
        titletext = 'SeaSoar tracks 1-6, layer: '+str(int(depthmin))+'-'+str(int(depthmax))+' m, '+str(frq)+' kHz'
        figbasename=figbasename0+'_SeaSoar'+str(int(depthmin))+'_'+str(int(depthmax))+'m'
    
    print gooddepth

    u2plot=u1[gooddepth,:]
    v2plot=v1[gooddepth,:]
    u2plot[abs(u2plot)>vlimit]=np.nan
    v2plot[abs(v2plot)>vlimit]=np.nan

    lon1[np.logical_and((lon1>-11.8),(lon1<-11.35))]=np.nan
      
    if doplot == 1:
        fig=plt.figure()
        ax = fig.add_subplot(111)
        llon,llat = m(lon1,lat1)
        Q=plt.quiver(llon,llat,u2plot,v2plot,\
                     units='width',width=0.002,color='k',alpha=0.9,label="HF-Radar",
                     headwidth=0,scale=5)

        # Put reference vector    
        k = plt.quiverkey(Q,0.75,0.15, 0.5, r'$0.5\, ms^{-1}$',labelpos='S',
                   fontproperties={'weight':'bold','size':'18'},color='k',zorder=3)
        plt.plot(llon,llat,'--',color="0.7",zorder=2,lw=0.3)
           
        m.drawcoastlines(ax=ax,zorder=3)
        m.fillcontinents(color='0.7', ax=ax,zorder=2)

        m.drawparallels(lat2plot, linewidth=0.5,labels=[1, 0, 0, 0] ,fontsize=14,ax=ax,zorder=1)
        m.drawmeridians(lon2plot, linewidth=0.5,labels=[0, 0, 0, 1] ,fontsize=14,ax=ax,zorder=1)
        
        plt.title(titletext,fontname='Times New Roman',fontsize=20)
        plt.savefig(figdir+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.show()
        plt.close()
  

