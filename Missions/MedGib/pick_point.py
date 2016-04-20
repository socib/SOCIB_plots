import test_pick_point
import numpy as np
import glob
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from matplotlib  import colors
from mpl_toolkits.basemap import Basemap
import medgib_functions

#points = test_pick_point.pic2data('/data_local/Publis/20150528_MedGib/figures/1-s2.0-S0079661113001080-gr1B.jpg')
points = test_pick_point.pic2data('/data_local/Publis/20150528_MedGib/figures/Medgib_SST_medium_domain_currents2_reversed.png')
np.savetxt('/home/ctroupin/Publis/20150528_MedGib/python/AtlanticJet2.dat', points)
