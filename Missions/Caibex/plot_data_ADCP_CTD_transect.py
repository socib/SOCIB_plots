#--------------------------------------------------------------------------
#
# plot_data_ADCP_CTD_transect.py
#
# plot the ADCP velocities and the CTD positions
#
# ctroupin, November 2013
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
basemap_resolution = 'l'
doplot = 1
frq = 75
vlimit = 2.0

# files and directories
ctddir  = '/home/ctroupin/ULPGC/CAIBEX_campaign/data/CTD/' 
datadir = '/home/ctroupin/ULPGC/CAIBEX_campaign/data/ADCP/' 
figdir = '/home/ctroupin/Publis/CAIBEX_CSR2012/figures/1/' 
ctdfile = 'CTDtransect_coord.dat' 
figdir = '/home/ctroupin/ULPGC/CAIBEX_campaign/data/ADCP/figures_python/CTD_transect/'
figtype='.png'

if not(os.path.exists(figdir)):
    os.makedirs(figdir)

figbasename = '5_adcp_'+str(frq)+'khz' 

# --------------------------------------------------------------------------------------

# region of interest
lonmin = -11 
lonmax = -9.5 
latmin = 30 
latmax = 31.3333 

# Spacing between x/y labels
dlon = 0.5
dlat = 0.5

# prepare projection
m = Basemap(projection='cyl',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                resolution=basemap_resolution)


# time for the CTD transect: 25-27 August 2009
timeref=dt.datestr2num('2009,1,1')
timestart1 = dt.datestr2num('2009,8,25')-timeref+1 
timeend1 = dt.datestr2num('2009,8,27')-timeref+1 

# depth for the plot
depthmin = 75
depthmax = 125

# --------------------------------------------------------------------------------------

# load CTD data
data=np.loadtxt(ctddir+ctdfile,usecols=(0,1,2))
lonCTD= data[:,0] 
latCTD = data[:,1] 
tempCTD = data[:,2] 

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

# --------------------------------------------------------------------------------------

# find the good depth

if (depthmin == depthmax):
    gooddepth = np.argmax(abs(zc-depthmin))
    u = u[gooddepth,:]
    v = v[gooddepth,:]
else:
    gooddepth = np.nonzero(np.logical_and(zc>=depthmin,zc<=depthmax))
    gooddepth = gooddepth[0]
    # average velocities over the depth range
    u = np.mean(u[gooddepth,:],axis=0)
    v = np.mean(v[gooddepth,:],axis=0)
 
# select good time
goodtime =  np.nonzero(np.logical_and(time>=timestart1,time<=timeend1)) 
goodtime = range(2714,3167)    
lon1 = lon[goodtime] 
lat1 = lat[goodtime] 
u1 = u[goodtime] 
v1 = v[goodtime] 


# Mask largest values (unrealistic)
goodval = np.nonzero(np.logical_and(abs(u1)<=vlimit,abs(v1)<=vlimit)) 
u1 = u1[goodval]
v1 = v1[goodval]
lon1 = lon1[goodval]
lat1 = lat1[goodval]

if depthmin == depthmax:
    titletext = 'CTD transect, depth: '+str(depthmin)+' m, '+str(frq)+' kHz'
    figbasename=figbasename+'_CTD_transect'+str(depthmin)+'m'
else:
    titletext = 'CTD transect, layer: '+str(depthmin)+'-'+str(depthmax)+' m, '+str(frq)+' kHz'
    figbasename=figbasename+'_CTD_transect'+str(depthmin)+'_'+str(depthmax)+'m'


# interpolate 
lat2interp = np.arange(30.2,31.2,0.01) 
lon2interp = -10.6*np.ones_like(lat2interp) 
#
u1interp = interpolate.griddata(np.transpose(np.array((lon1,lat1))),u1,
                                np.transpose(np.array((lon2interp,lat2interp))),method='linear')
v1interp = interpolate.griddata(np.transpose(np.array((lon1,lat1))),v1,
                                np.transpose(np.array((lon2interp,lat2interp))),method='linear')


##newxtick = [-11 lon2interp(1) -10] 
##newytick = 30:0.5:31 
##
# compute Coriolis frequency
f = 2*7.2921e-5*np.sin(np.deg2rad(30.5))

if doplot == 1:
    fig=plt.figure()
    ax = fig.add_subplot(111)
    lon2interp,lat2interp = m(lon2interp,lat2interp)
    Q=plt.quiver(lon2interp,lat2interp,u1interp,v1interp,\
                 units='width',width=0.003,color='k',alpha=0.9,label="HF-Radar")

    # Put reference vector    
    k = plt.quiverkey(Q,0.5,0.2, 0.2, r'$0.2\, ms^{-1}$',labelpos='S',
               fontproperties={'weight':'bold','size':'18'},color='k')

    # Add CTD points
    m.plot(lonCTD,latCTD,'ko',ms=7,color='.6') 
    
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='black', ax=ax)

    m.drawparallels(np.arange(latmin,latmax,dlat), linewidth=0.5,
                        labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16)
    m.drawmeridians(np.arange(lonmin,lonmax,dlon), linewidth=0.5,
                        labels=[0, 0, 0, 1], fontname='Times New Roman',fontsize=16)
    
    plt.title(titletext,fontname='Times New Roman',fontsize=20)
    plt.savefig(figdir+figbasename+figtype, dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)
    #plt.show()
    plt.close()
  

##    delta_u = u1interp(2:end) - u1interp(1:end-1) 
##    delta_y = 6371000*deg2rad(distance(lat2interp(2:end),...
##    lon2interp(2:end),lat2interp(1:end-1),lon2interp(1:end-1))) 
##    shear = -delta_u./delta_y 
##    relative_shear = shear/f 
##    
##      
##   ##    
##    [llon,llat] = m_ll2xy(0.5*(lon2interp2(1:end-1)+lon2interp2(2:end)),...
##                          0.5*(lat2interp(1:end-1)+lat2interp(2:end))) 
##    scatter(llon,llat,10,relative_shear,'filled','marker','s') 
##    scatter(llon-0.0002,llat,10,relative_shear,'filled','marker','s') 
##    scatter(llon-0.0004,llat,10,relative_shear,'filled','marker','s') 
##    scatter(llon-0.0006,llat,10,relative_shear,'filled','marker','s') 
##    hc = colorbar 
##    set(hc,'pos',[0.5121 0.2452 0.03259 0.5381]) 
##    set(get(hc,'title'),'string','Shear/f','fontsize',16) 
##    colormap(cmap) 
##    
##    
##    
##    m_text(-10.73,latCTD(end),'T_{13}','fontsize',14,...
##        'fontweight','bold','verticalalign','c') 
##    m_text(-10.73,latCTD(1),'T_{1}','fontsize',14,...
##        'fontweight','bold','verticalalign','c') 
##    m_text(-10.8,latCTD(5),'T_{5}','fontsize',14,...
##        'fontweight','bold','verticalalign','c') 
##
##
##
### compute max velocity
##unorm = sqrt(u1.*u1+v1.*v1) 
##min(unorm)
##max(unorm)
