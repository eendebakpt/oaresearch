# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 11:47:57 2012

@author: eendebakpt
"""


# %% Load necessary packages """
from __future__ import print_function

import os
import numpy as np
import matplotlib.pyplot as plt
oadir = '/home/eendebakpt/misc/oa/oacode/'
import oapackage
from oapackage import *


def tickfontsize(fontsize=14, ax=None):
    if ax is None:
        ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    plt.draw()


def nonregularProperties(sols):
    D = [A.Defficiency() for A in sols]

    nn = len(sols)
    A3 = np.zeros(nn, )
    A4 = np.zeros(nn, )

    for i, al in enumerate(sols):
        g = al.GWLP()
        A3[i] = g[3]
        A4[i] = g[4]
    return D, A3, A4


# %% Load data

if 1:
    # for paper
    arrayclass = oapackage.arraydata_t(2, 36, 2, 7)
    #arrayclass=oapackage.arraydata_t(2, 40, 2, 7)

    basedir = '/home/eendebakpt/misc/homepage/oapage/tpages/'
    basedir = '/home/eendebakpt/misc/oapage2/tpages/'
    oafile = 'class-%s-t%dselection.oa' % (arrayclass.idstr(), arrayclass.strength)
    afile = os.path.join(basedir, 'classdata-special-%s-t%d' % (arrayclass.idstr(), arrayclass.strength), oafile)
else:
    arrayclass = oapackage.arraydata_t(2, 36, 2, 5)   # for paper
    arrayclass = oapackage.arraydata_t(2, 40, 2, 5)  # for paper
    #arrayclass=oapackage.arraydata_t(2, 56, 2, 5)

    basedir = '/home/eendebakpt/misc/homepage/oapage/tpages/'
    oafile = 'class-%s-t%dselection.oa' % (arrayclass.idstr(), arrayclass.strength)
    afile = os.path.join(basedir, 'abdata-%s-t%d' % (arrayclass.idstr(), arrayclass.strength), oafile)

sols = oalib.readarrayfile(afile)

print('read %d arrays' % len(sols))

# %% Calculate properties

D, A3, A4 = nonregularProperties(sols)


# %%

def formatFigure(fig=None, paperfig=True):
    if fig is None:
        fig = plt.gcf()
    else:
        plt.figure(fig)
    if paperfig:
        plt.title('')
# %% Show


lstr = arrayclass.latexstr()
tstr = 'Selection of %d arrays in $%s$' % (len(sols), lstr)

plt.figure(10)
plt.clf()
plt.plot(A3, D, '.b', markersize=12)
plt.xlabel('$A_3$', fontsize=14)
plt.ylabel('D-efficiency', fontsize=14)
plt.title(tstr, fontsize=17)

plt.figure(11)
plt.clf()
plt.plot(A3 + A4, D, '.b', markersize=12)
plt.xlabel('$A_3+A_4$', fontsize=14)
plt.ylabel('D-efficiency', fontsize=14)
plt.title(tstr, fontsize=17)


plt.figure(12)
plt.clf()
plt.plot(A3, A4, '.b', markersize=12)
plt.xlabel('$A_3$', fontsize=14)
plt.ylabel('$A_4$', fontsize=14)
plt.title(tstr, fontsize=17)


rcolor = [.8, 0, 0]
bcolor = [0, .3, 1]

fw = 'medium'
fw = 'normal'
# fw='semibold'
fig = plt.figure(20)
plt.clf()
plt.plot(D, A3, '.', color=bcolor, markersize=12, label='$A_3$')
plt.plot(D, A3 + A4, '.', color=rcolor, markersize=12, label='$A_3+A_4$')
plt.xlabel('D-efficiency', fontsize=22, fontweight=fw)
plt.ylabel('$A_3$, $A_3+A_4$', fontsize=22, fontweight=fw)
plt.title(tstr, fontsize=22)
ax = plt.gca()
tickfontsize(16)
legendh = plt.legend(numpoints=1, fontsize=18)

oapackage.niceplot(ax, fig=fig, legend=legendh)
formatFigure()

tilefigs([10, 11, 12, 20], [2, 2])


# %%

plt.title('')

import tempfile
picturedir = tempfile.mkdtemp(prefix='pictures-A3-A4-D')
idstr = arrayclass.idstr().replace('.', '-d-') + '-t%d' % arrayclass.strength
plt.figure(10)
plt.savefig(os.path.join(picturedir, 'A3-D-%s.png' % idstr))
plt.figure(11)
plt.savefig(os.path.join(picturedir, 'A34-D-%s.png' % idstr))
plt.figure(20)
plt.savefig(os.path.join(picturedir, 'A3-A34-%s.png' % idstr))
print('written figures to %s' % picturedir)

# %%
