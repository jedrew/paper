from __future__ import print_function, division
from astropy.io import fits, ascii
import numpy as np
import matplotlib as plt
from matplotlib.colors import LogNorm
#plt.rc('text', usetex=True)                                                    
#plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rcParams['text.latex.unicode']=True

from pylab import *

WIDTH = 441.01
FACTOR = 1.25
fig_width_pt  = WIDTH * FACTOR
inches_per_pt = 1.0 / 72.27
golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio  # figure height in inches
fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list

def generate(fn):
    basefn = fn.split('/')[-1].split('.')[0]
    data = fits.getdata(fn, 1)
    mask = np.where(data['r']-data['r-i']+data['Ai_SS']<19)
#    fig = figure(figsize=fig_dims)
    fig = figure()
    ax = fig.add_subplot(111)
    #imshow(dmap_i[btoj(bmin):btoj(bmax),ltoi(lmin):ltoi(lmax)], origin='lower', extent=[lmin,lmax,bmin,bmax], vmin=0, vmax=vmax, cmap='YlGnBu_r')
    #ax.xaxis.set_major_formatter(x_formatter)
    #ax.yaxis.set_major_formatter(y_formatter)
    ylim(25, 8)
    xlim(-0.5, 2.5)
    ax.scatter(data['r-i'], data['i'], c='lightgray', s=2, edgecolor='none')
    ax.scatter(data['r-i'][mask], data['i'][mask], c='r', s=2, edgecolor='none')
    ylabel('$r^{\prime}$')
    xlabel('$r^{\prime}-i^{\prime}$')
    subplots_adjust(bottom=0.175, top=0.925, right=0.955, left=0.125)
    savefig('../figures/cmd-{0}.png'.format(basefn), dpi=200)

generate('../data/1413061465.506936.resu.ext.fits')
generate('../data/1414018314.240538.resu.ext.fits')
