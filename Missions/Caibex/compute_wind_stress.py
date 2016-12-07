import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
import matplotlib.dates as dt
from matplotlib  import colors
from scipy import interpolate
from mpl_toolkits.basemap import Basemap

# define constants

g = 9.8        # acceleration due to gravity [m/s^2]
sigmaSB     = 5.6697e-8  # Stefan-Boltzmann constant [W/m^2/K^4]
eps_air     = 0.62197    # molecular weight ratio (water/air)
gas_const_R = 287.04     # gas constant for dry air [J/kg/K]
CtoK        = 273.16     # conversion factor for [C] to [K]

kappa          = 0.4     # von Karman's constant
Charnock_alpha = 0.011   # Charnock constant (for determining roughness length
R_roughness   = 0.11     # limiting roughness Reynolds # for aerodynamically

                      
# ------ defaults suitable for boundary-layer studies
cp            = 1004.7    # heat capacity of air [J/kg/K]
rho_air       = 1.22      # air density (when required as constant) [kg/m^2]
Ta_default    = 10        # default air temperature [C]


Ta=Ta_default

def [cd,u10]=cdntc(sp,z,Ta)

# iteration endpoint
tol=.00001

# Compute air viscosity from temperature
visc=vis = 1.326e-5*(1 + 6.542e-3*Ta + 8.301e-6*Ta**2 - 4.84e-9*Ta**3) 

### Allocate
ustaro=np.zeros_like(sp) 
ustarn=.036*sp 

# iterate to find z0 and ustar
ii=abs(ustarn-ustaro)>tol 
##while any(ii(:)),
##0043   ustaro=ustarn 
##0044   z0=Charnock_alpha.*ustaro.^2./g + R_roughness*visc./ustaro 
##0045   ustarn=sp.*(kappa./log(z./z0)) 
##0046   ii=abs(ustarn-ustaro)>tol 
##0047 end
##0048 
##0049 sqrcd=kappa./log((10)./z0) 
##0050 cd=sqrcd.^2 
##0051 
##0052 u10=ustarn./sqrcd 
##
##
##
##
##windstress = rhoa*(cd.*u10.^2)
