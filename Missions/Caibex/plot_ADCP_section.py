#--------------------------------------------------------------------------
# plot_ADCP_section.py
#
# plot velocity
#
# ctroupin, January 2014
#--------------------------------------------------------------------------

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
from scipy import interpolate
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import scipy.io

ADCPdir='/home/ctroupin/ULPGC/CAIBEX_campaign/data/ADCP/'
ADCPcoord1='CG75_xy.mat'
ADCPcoord2='CG150_xy.mat'
ADCPvel1='CG75_uv.mat'
ADCPvel2='CG150_uv.mat'

cmap_redblue=plt.cm.RdBu

goodtime = range(2715,3168)

velmin,velmax=-0.5,0.5
levels2plot3=np.arange(velmin,velmax,0.05)
normvel=colors.Normalize(velmin,velmax)

# Load data
mat = scipy.io.loadmat(ADCPdir+ADCPcoord1)
xyt = mat['xyt']
depth= np.squeeze(mat['zc'])
lonADCP1 = -360.+xyt[0,goodtime]
latADCP1 = xyt[1,goodtime]
timeADCP1 = xyt[2,goodtime]

mat = scipy.io.loadmat(ADCPdir+ADCPcoord2)
xyt = mat['xyt']
depth2= np.squeeze(mat['zc'])
lonADCP2 = -360.+xyt[0,goodtime]
latADCP2 = xyt[1,goodtime]
timeADCP2 = xyt[2,goodtime]

mat = scipy.io.loadmat(ADCPdir+ADCPvel1)
uv = mat['uv']
uADCP1 = uv[:,0::2]
vADCP1 = uv[:,1::2]
uADCP1,vADCP1 = uADCP1[:,goodtime],vADCP1[:,goodtime]
mat = scipy.io.loadmat(ADCPdir+ADCPvel2)

uv = mat['uv']
uADCP2 = uv[:,0::2]
vADCP2 = uv[:,1::2]
uADCP2,vADCP2 = uADCP2[:,goodtime],vADCP2[:,goodtime]

# Mask
uADCP1=np.ma.masked_where(np.isnan(uADCP1),uADCP1)
vADCP1=np.ma.masked_where(np.isnan(vADCP1),vADCP1)
uADCP2=np.ma.masked_where(np.isnan(uADCP2),uADCP2)
vADCP2=np.ma.masked_where(np.isnan(vADCP2),vADCP2)

# plot
fig=plt.figure(figsize=(10,5))
ax=fig.add_subplot(121)
plt.contourf(latADCP1,-depth,uADCP1,levels2plot3,norm=normvel,cmap=cmap_redblue)
ax=fig.add_subplot(122)
plt.contourf(latADCP2,-depth2,uADCP2,levels2plot3,norm=normvel,cmap=cmap_redblue)
plt.show()
plt.close()
