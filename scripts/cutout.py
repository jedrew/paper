from __future__ import print_function, division
from astropy.io import fits, ascii
import numpy as np
import matplotlib as plt
from matplotlib.colors import LogNorm
plt.rc('text', usetex=True)                                                    
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams['text.latex.unicode']=True

from pylab import *

def btoj(b):
    return (b + 6) * 60
def ltoi(l):
    return (l - 25) * 60 

def fmtb(x, x2):
    if x<=0:
        return u"{0:.0f}\u00b0".format(x)
    else:
        return u"+{0:.0f}\u00b0".format(x)

x_formatter = plt.FormatStrFormatter(u"%.0f\u00b0") 
y_formatter = FuncFormatter(fmtb)

WIDTH = 441.01
FACTOR = 1.25
fig_width_pt  = WIDTH * FACTOR
inches_per_pt = 1.0 / 72.27
golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio  # figure height in inches
fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list

dmap_i = fits.getdata("../data/dmap-i-1-cube.fits", 3)

def generate(lmin, lmax, bmin, bmax, vmax):
    fig = figure()
    ax = fig.add_subplot(111)
    imshow(dmap_i[btoj(bmin):btoj(bmax),ltoi(lmin):ltoi(lmax)], origin='lower', extent=[lmin,lmax,bmin,bmax], vmin=0, vmax=vmax, cmap='YlGnBu_r')
    ax.xaxis.set_major_formatter(x_formatter)
    ax.yaxis.set_major_formatter(y_formatter)
    xlabel('Galactic longitude')
    ylabel('Galactic latitude')
    colorbar(label='Sources per sq. arcmin.', extend='max')
    savefig('../figures/dmap-cutout-{0}-{1}-{2}-{3}.pdf'.format(lmin, lmax, bmin, bmax))

generate(30, 40, -5, 5, 75)
generate(80, 90, -5, 5, 30)
generate(170,180, -5, 5, 15)

