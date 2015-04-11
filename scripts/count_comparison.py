from __future__ import print_function, division
from astropy.io import fits, ascii
import iphas2gaia
import numpy as np
import matplotlib as plt
plt.rc('text', usetex=True)                                                    
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams['text.latex.unicode']=True

from pylab import *
import sys
sys.path.append('/local/home/hfarnhill/PycharmProjects/extinction/')
import sightline
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def get_values(e, e_err):
    err = 1.0 / np.sum(np.power(e_err, -2))
    av_e = np.sum(np.array(e)/np.power(np.array(e_err),2)) * err
    return av_e, err

def get_l(i):
    return (l - 25) / (12.0/60.0)

def get_b(j):
    return (b +6) / (12.0/60.0)

def fmtb(x, x2):
    if x<=0:
        return u"{0:.0f}\u00b0".format(x)
    else:
        return u"+{0:.0f}\u00b0".format(x)
y_formatter = FuncFormatter(fmtb)

counts = fits.getdata(sys.argv[1])
l0 = float(sys.argv[2])
ymax=int(sys.argv[3])

fig = figure(figsize=(8,5))

axb = fig.add_subplot(111)
axb.spines['top'].set_color('none')
axb.spines['bottom'].set_color('none')
axb.spines['left'].set_color('none')
axb.spines['right'].set_color('none')
axb.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
axb.set_xlabel('Galactic latitude')
axb.set_ylabel(r'Source density (per sq. arcmin.)', labelpad=10)

types = ['m', 's', 'l', 'SS', 'i']
cs = ['r', 'b', 'g', 'orange', 'k']

types = ['i', 'SS', 'l', 's', 'm']
cs = ['k', 'orange', 'g', 'b', 'r']

types = ['i', 'l', 'SS',  's', 'm']
cs = ['k', 'g', 'orange', 'b', 'r']

types = ['l', 'SS',  's', 'm']
cs = ['g', 'orange', 'r', 'b']

axes = {}

for i in range(1):
    lmin = l0 + i
    if 85<lmin<91:
        fac = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        fac2 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        majx = 10
        types = ['l', 'SS',  's', 'm']
        cs = ['g', 'orange', 'r', 'b']
        addx = 0
        bottom = 0.085
    elif 170<lmin<180:
        fac = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        fac2 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        majx = 2
        types = ['l', 'SS',  's' ]
        cs = ['g', 'orange', 'r' ]
        addx = 0.2
        bottom = 0.12
    elif 28<lmin<35:
#        fac = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
#        fac = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        fac = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        fac2 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        majx = 20
        types = ['l', 'SS', 'm']
        cs = ['g', 'orange', 'b']
        addx = 0.2
        bottom = 0.12

    print(fac)

    #ax = fig.add_subplot(2, 1, 1 + i)
    fig = figure(figsize=(8,1.25*len(types)+addx))

    axb = fig.add_subplot(111)
    axb.spines['top'].set_color('none')
    axb.spines['bottom'].set_color('none')
    axb.spines['left'].set_color('none')
    axb.spines['right'].set_color('none')
    axb.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    axb.set_xlabel('Galactic latitude')
    axb.set_ylabel(r'Source density (per sq. arcmin.)', labelpad=10)
    for j, t in enumerate(reversed(types)):
        axes[j] = fig.add_subplot(len(types), 1, 1 + j)
        if l0>90 and t=='m':
            continue
        if l0<65 and t=='s':
            continue
        mask = counts['counts_i_{1}'.format(t, lmin)] > 0
        print(len(counts))
        print(len(mask))
        axes[j].plot(counts['b'][mask], counts['counts_i_{1}'.format(t, lmin)][mask]/144.0, c='k')
        axes[j].scatter(counts['b'][mask], counts['counts_i_{1}'.format(t, lmin)][mask]/144.0, edgecolor='none', s=10, c='k')
        mask = counts['counts_{0}_{1}'.format(t, lmin)] > 0
        temp = counts[mask]
        if t!='l':
            axes[j].plot(counts['b'][mask], counts['counts_{0}_{1}'.format(t, lmin)][mask]/(144.0*fac[mask]), c=cs[len(types)-1-j])
            axes[j].scatter(counts['b'][mask], counts['counts_{0}_{1}'.format(t, lmin)][mask]/(144.0*fac[mask]), edgecolor='none', s=10, c=cs[len(types)-1-j])
        else:
            axes[j].plot(counts['b'][mask], counts['counts_{0}_{1}'.format(t, lmin)][mask]/(144.0*fac2[mask]), c=cs[len(types)-1-j])
            axes[j].scatter(counts['b'][mask], counts['counts_{0}_{1}'.format(t, lmin)][mask]/(144.0*fac2[mask]), edgecolor='none', s=10, c=cs[len(types)-1-j])
        axes[j].xaxis.set_major_formatter(y_formatter)
        axes[j].yaxis.set_major_locator(MultipleLocator(majx))
        axes[j].xaxis.set_minor_locator(MultipleLocator(0.5))
        axes[j].xaxis.set_major_locator(MultipleLocator(1))
        if j!=len(types)-1:
            axes[j].set_xticklabels([])
        ylims = axes[j].get_ylim()
        ylim(0, ylims[1])
        ylim(0, ymax)
        xlim(-5.5,5.5)


maxes = []
mins = []
for i in range(4):
    maxes.append(np.max(axes[j].get_ylim()))
    mins.append(np.min(axes[j].get_ylim()))
absmin = np.min(mins)
absmax = np.max(maxes)
for i in range(4):
    axes[j].set_ylim(absmin, absmax)
    axes[j].update({'ylim':(absmin, absmax)})
    print(axes[j].get_ylim())

draw()
#show()

subplots_adjust(top=0.95, bottom=bottom, left=0.08, right=0.975)
fig.savefig('count_comparison_{0}.pdf'.format(int(l0)))
fig.savefig('count_comparison_{0}.png'.format(int(l0)))
#show()


