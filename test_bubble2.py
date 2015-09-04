from pylab import *
import numpy as np
from matplotlib.patches import Ellipse

doplot=1
# reading the data from a csv file
datafile = '/home/ctroupin/SOCIB/Python/Plot_scales/timescale3.txt'
labels = np.genfromtxt(datafile, usecols=(0, ), delimiter=',', dtype=str)
rdata = np.genfromtxt(datafile, usecols=(1, 2, 3, 4), delimiter=',')
x=rdata[:, 0]
y=rdata[:, 1]
ellwidth=rdata[:, 2]
ellheight=rdata[:, 3]

if doplot:
    fig=plt.figure(num=None, figsize=(10,6), facecolor='w')
    ax = fig.add_subplot(111)
    # making the scatter plot

    for j in xrange(0, 4):
        ellipse = Ellipse(xy=(x[j], y[j]), width=2.*ellwidth[j], height=2.*ellheight[j], 
                        edgecolor='r', fc='None', lw=2)
        ax.add_patch(ellipse)

    axis([0,11,0,11])
    xlabel('Horizontal spatial scales')
    ylabel('Time scales')
    ax.set_xticks(np.arange(0, 10))
    ax.set_yticks(np.arange(0, 10))
    plt.show()
    plt.close()
