"""
Created on Fri Nov 16 11:24:39 2012

@author: eendebakpt
"""


from researchOA import *
import researchOA
import oalib
import operator
import matplotlib.pyplot as plt
import time
import numpy as np
import os
import sys
oadir = '/home/eendebakpt/misc/oa/oacode/'
xdir = '/home/eendebakpt/misc/oa/dof'

""" Load necessary packages """
sys.path.append(os.path.join(oadir, 'pythondevelop'))
sys.path.append(os.path.join(oadir, 'oalib'))

os.chdir(xdir)
t = 2
N = 16
ncols = 8
kstart = 6
aidx = 16
N = 16
ncols = 8
kstart = 7
aidx = 27  # forward order
N = 16
ncols = 8
kstart = 6
aidx = 3  # forward order
N = 16
ncols = 8
kstart = 7
aidx = 4  # forward order
# N=32; t=3; ncols=9;kstart=9; aidx=7 # forward order
# N=8; ncols=8;kstart=4; aidx=0 # reverse order
# N=12; ncols=8;kstart=5; aidx=0 # ??
adata = oalib.arraydata_t(oalib.intVector([2] * ncols), N, t, ncols)
adata0 = oalib.arraydata_t(adata, kstart)
adata0.writeConfigFile('oaconfig%d.txt' % kstart)
afile0, nsols = extendInitial(
    xdir, adata, kstart, verbose=2, cache=0, cmdlog=None)


sols = oalib.readarrayfile(afile0)
al = sols[aidx]
aldof = oalib.reduceDOPform(al)
al.showarray()
aldof.showarray()

X = al.getarray()
X.T.flatten()


# math.factorial(9)

# binom(9,5)*math.factorial(5)

def showVector(v):
    for ii in range(0, v.size()):
        print(('%s,' % v[0]), end='')
    print('')


print('')
print('')


ncols = 4
ad = oalib.arraydata_t(2, 8, 2, ncols)
oaextendoptions = oalib.OAextend(ad)
al = oalib.array_link(8, 2, 0)
al.create_root(ad)
al.showarray()
alist = oalib.arraylist_t()
alist.push_back(al)
alist2 = oalib.extend_arraylist(alist, ad, oaextendoptions)
for al in alist2:
    al.showarray()
