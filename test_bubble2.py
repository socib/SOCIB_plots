from pylab import *
import numpy as np
from matplotlib.patches import Ellipse

figdir = '/home/ctroupin/Presentations/figures4presentations/'
figname = 'scales_python'

doplot=1
# reading the data from a csv file
datafile = './timescale3.txt'
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

    for j in xrange(0, 5):
        ellipse = Ellipse(xy=(x[j], y[j]), width=1.*ellwidth[j], height=1.*ellheight[j],
                        edgecolor='.8', fc='None', lw=2)
        ax.add_patch(ellipse)
        plt.text(x[j], y[j], labels[j], horizontalalignment='center', verticalalignment='center')

    axis([0,11,0,11])
    xlabel('Horizontal spatial scales')
    ylabel('Time\n scales', rotation=0, ha='right')
    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 8))
    plt.xlim(0,7.5)
    plt.ylim(0,8.5)
    ax.set_xticklabels(('10 m','100 m','1 km','10 km','100 km','1000 km', '10000 km', '100000 km'))
    ax.set_yticklabels(('1 hour','1 day','1 week','1 month','1 year','10 years', '100 years', '1000 years'))
    plt.show()
    plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close()
