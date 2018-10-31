# -*- coding: utf-8 -*-
"""

Example script for calculating conference designs
    
Pieter Eendebak <pieter.eendebak@gmail.com>
"""

#%% Load necessary packages
import os
import platform
import sys
import numpy as np
from os.path import join
import itertools

import oapackage
from oapackage import oahelper

import oapackage.graphtools

oadir = os.path.join(os.path.expanduser('~'), 'misc/oa/oacode/')
sys.path.append(os.path.join(oadir, 'pythondevelop'))
from researchOA import job, createJobScript

#%% Setup directories
basedir=os.path.expanduser('~')
resultsdir = join(basedir, 'oatmp')
outputdir = oapackage.mkdirc(
    os.path.join(os.path.expanduser('~'), 'oatmp', 'conf2'))

scriptdir = join(os.path.expanduser('~'), 'confjobs')

full_generation = False
verbose=1


def reportScriptFile(scriptfile, verbose=1):
    nlines=0
    ncommands=0
    with open(scriptfile, 'rt') as fid:
        for line in fid:
            if verbose>=2:
                print(line.strip())
            nlines=nlines+1
            if not line.startswith('#'):
                ncommands=ncommands+1
    print('generated script file: %s with %d commands'  % (scriptfile, nlines))
    
        
#%%

def generateConference(N, kmax=None, verbose=1, diagc=False, nmax=None, selectmethod='random', tag='cdesign', outputdir=None):
    """ General function to compute conference matrices 

    Arguments:
        N : integer
            number of rows in the array
        kmax : integer
            maximum number of columns to compute
        verbose : integer
            output level
        diagc : boolean
            the default value is False. If True, then only the diagonal
            matrices will be computed (e.g. all zeros are on the diagonal)

    """
    if kmax is None:
        kmax = N
    ctype = oapackage.conference_t(N, N, 0)

    if diagc:
        ctype.ctype = oapackage.conference_t.CONFERENCE_DIAGONAL
        tag += '-diagonal'
    if not nmax is None:
        tag += '-r'

    al = ctype.create_root()

    ll = oapackage.arraylist_t()
    ll.push_back(al)
    LL = [[]] * (kmax)
    LL[1] = ll
    print('generateConference: start: %s' % ctype)
    if outputdir is not None:
        _ = oapackage.writearrayfile(
            join(outputdir, 'cdesign-%d-%d.oa' % (N, 2)), LL[1], oapackage.ATEXT, N, 2)

    for extcol in range(2, kmax):
        if verbose:
            print('generateConference: extcol %d: %d designs' % (extcol, len(LL[extcol - 1])) )
            sys.stdout.flush()
        LL[extcol] = oapackage.extend_conference(
            LL[extcol - 1], ctype, verbose=verbose >= 2)

        LL[extcol] = oapackage.selectConferenceIsomorpismClasses(LL[extcol], 1)

        LL[extcol] = oapackage.sortLMC0(LL[extcol])

        if nmax is not None:
            na = min(nmax, len(LL[extcol]))
            if na > 0:
                if selectmethod == 'random':
                    idx = np.random.choice(len(LL[extcol]), na, replace=False)
                    LL[extcol] = [LL[extcol][i] for i in idx]
                elif selectmethod == 'first':
                    LL[extcol] = [LL[extcol][i] for i in range(na)]
                else:
                    # mixed
                    raise Exception('not implemented')
        afmode = oapackage.ATEXT
        if (len(LL[extcol]) > 1000):
            afmode = oapackage.ABINARY
        if outputdir is not None:
            _ = oapackage.writearrayfile(join(
                outputdir, '%s-%d-%d.oa' % (tag, N, extcol + 1)), LL[extcol], afmode, N, extcol + 1)

    ll = [len(l) for l in LL]
    if verbose:
        print('generated sequence: %s' % ll)
    return LL


LL = generateConference(12, outputdir=outputdir)


#%% Test maxz values


def showMaxZ(LL):
    N = LL[3][0].n_rows

    for ii, L in enumerate(LL):
        k = ii + 1
        s = [oapackage.maxz(al) for al in L]
        mm, _ = np.histogram(s, range(N + 1))
        print('%d cols: maxz seq %s' % (k, list(mm)) )

showMaxZ(LL)


#%%
if not full_generation:
    # only do small cases
    for NN in range(4, 12, 2):
        _ = generateConference(N=NN)

else:
    for NN in range(4, 18, 2):
        _ = generateConference(N=NN)
    for NN in range(20, 26, 2):
        _ = generateConference(N=NN, kmax=3)
    # big cases (takes a longer time)
    LL = generateConference(N=18)
    LL = generateConference(N=20, kmax=20)
    LL = generateConference(N=22, kmax=7)
    _ = generateConference(N=24, kmax=6)
    LL = generateConference(N=26, kmax=6)
    LL = generateConference(N=28, kmax=5)
    LL = generateConference(N=30, kmax=4)

    LL = generateConference(22, diagc=True)
    LL = generateConference(24, diagc=True)

    LL = generateConference(26, diagc=True)
    LL = generateConference(28, kmax=9, diagc=True)  # started...
    #LL=generateConference(30, kmax=7, diagc=True)

    LL = generateConference(28, diagc=True, nmax=100)

    cmd = './oaconference  -N 22  -o cdesign'
    print('do manually: %s' % cmd)
    cmd = './oaconference  -N 22 -k 7 -o cdesign'
    print('do manually: %s (exceeds memory for next column)' % cmd)


#%%
scriptfile = os.path.join(scriptdir, 'jobfile-small.sh')
with open(scriptfile, 'wt') as fid:
    for N in range(2, 18, 2):
        cmd = 'oaconference -N %d --ctype 0 --itype 2  --j1zero 0 --j3zero 0 -o cdesign' % N  # conference matrix
        _ = fid.write('%s\n' % cmd)

        # double conference design
        cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 1 --j3zero 1 -o dconferencej1j3' % N
        _ = fid.write('%s\n' % cmd)
        cmd = 'oaconference -N %d --ctype 2 --itype 2  --j1zero 0 --j3zero 0 -o weighing2' % N
        _ = fid.write('%s\n' % cmd)
        #cmd='oaconference -N %d --ctype 2 --itype 3  --j1zero 0 --j3zero 0 -o dconference2' % N
        # print(cmd)
        #_=fid.write('%s\n'  % cmd)

        cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 1 --j3zero 0 -o dconferencej1' % N
        _ = fid.write('%s\n' % cmd)
        
        # J3=0, J1 unrestricted
        cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 0 --j3zero 1  -o dconferencej3' % N
        _ = fid.write('%s\n' % cmd)

    for (N, kmax) in [(20, 4), (18, 4), (22, 4), (24, 2)]:
        cmd = 'oaconference -N %d -k %d --ctype 2 --itype 2  --j1zero 0 --j3zero 0 -o weighing2' % (N, kmax)
        _ = fid.write('%s\n' % cmd)

    for N in [18, 20]:
        cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 1 --j3zero 0 -o dconferencej1' % N
        _ = fid.write('%s\n' % cmd)
    for N in [26]:
        cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 1 --j3zero 0 -k 5 -o dconferencej1' % N
        _ = fid.write('%s\n' % cmd)

    for N in [18, 20]:
        cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 0 --j3zero 1  -o dconferencej3' % N
        _ = fid.write('%s\n' % cmd)

    for N in [30, 32, 34]:
        cmd = 'oaconference -N %d --ctype 0 --itype 2  --j1zero 0 --j3zero 0 -o cdesign -k 5' % N  # conference matrix
        _ = fid.write('%s\n' % cmd)
        cmd = 'oaconference -N %d --ctype 1 --itype 2  --j1zero 0 --j3zero 0 -o cdesign-diagonal -k 7' % N  # conference matrix
        _ = fid.write('%s\n' % cmd)

    for N in [34]:
        cmd = 'oaconference -N %d --ctype 1 --itype 2  --j1zero 0 --j3zero 0 --select 2 -o cdesign-diagonal-r -k 7' % N  # conference matrix
        _ = fid.write('%s\n' % cmd)

    for N in [40]:
        cmd = 'oaconference -N %d --ctype 0 --itype 2  --j1zero 0 --j3zero 0 -o cdesign -k 5' % N  # conference matrix
        _ = fid.write('%s\n' % cmd)

    # random extension cdesigns
    for N in [30]:
        cmd = 'oaconference -N %d --ctype 1 --itype 2  --j1zero 0 --j3zero 0 -o cdesign-diagonal-r --select 2 --nmax 2000' % N  # conference matrix
        _ = fid.write('%s\n' % cmd)

reportScriptFile(scriptfile, verbose=1)

#%% Double conference matrices

if 1:
    scriptfile =os.path.join(scriptdir, 'jobfile.sh')
    with open(scriptfile, 'wt') as fid:
        for n in range(2, 12, 2):
            cmd = 'oaconference -N %d --ctype 2 --itype 3  --j1zero 0  --j3zero 0 -o dconference2' % n
            _ = fid.write('%s\n' % cmd)
        for (n, kmax) in [(14, 4), (12, 4), (16, 3), (18,3), (20,3)]:
            cmd = 'oaconference -N %d -k %d --ctype 2 --itype 3  --j1zero 0 --j3zero 0 -o dconference2' % (n, kmax)
            _ = fid.write('%s\n' % cmd)

        for n in range(18, 38, 2):
            cmd = 'oaconference -N %d --ctype 2 --itype 3 --j1zero 1 --j3zero 1  -o dconferencej1j3' % n
            _ = fid.write('%s\n' % cmd)

        for n in range(38, 45, 2):
            cmd = 'oaconference -N %d --ctype 2 --itype 3 --j1zero 1 --j3zero 1  -o dconferencej1j3' % n
            _ = fid.write('%s\n' % cmd)

        ww = [(46, 7), (48, 7), (50, 7), (52, 6)]
        ww += [(56, 6), (58, 6), (60, 5), (62, 24)]
        ww += [(n, 5) for n in range(64, 70, 2)] + [(70, 20)] + [(n, 5) for n in range(72, 82, 2)] + [(78, 20)]
        for n, kmax in ww:
            cmd = 'oaconference -N %d -k %d --ctype 2 --itype 3 --j1zero 1 --j3zero 1  -o dconferencej1j3' % (n, kmax)
            _ = fid.write('%s\n' % cmd)
    reportScriptFile(scriptfile)
    
    cmd = 'oaconference -N 48 --ctype 2 --itype 3  --j1zero 1 --j3zero 1 -o dconferencej1j3-r --nmax 400'
    print(cmd)

    cmd = 'oaconference -N 20 --ctype 0 --itype 2 --j1zero 1 --j3zero 0 -o cdesign-r --nmax 40'
    print(cmd)


    # partial enumeration double conf
    # nmaxs are chosen such that a single-core job can run within a day
    NN = [44, 48, 50, 52, 54, 56, 58]
    nmaxs = [200, 60000, 18000, 18000, 6000, 6000, 4000]
    jj = []
    for i, N in enumerate(NN):
        nmax = nmaxs[i]
        cmd = 'oaconference -N %d --ctype 2 --itype 3 --j1zero 1 --j3zero 1 -o dconferencej1j3-r --nmax %d' % (N, nmax)
        if verbose>=2:
            print(cmd)
        j = job(cmd, shorttag='N%d' % N)
        jj += [j]
    NN = [60, 62, 64, 66]
    nmaxs = [2000, 2000, 2000, 200]
    for i, N in enumerate(NN):
        nmax = nmaxs[i]
        cmd = 'oaconference -N %d --ctype 2 --itype 3 --j1zero 1 --j3zero 1 -o dconferencej1j3-r --nmax %d' % (N, nmax)
        if verbose>=2:
            print(cmd)
        j = job(cmd, shorttag='N%d' % N)
        jj += [j]

    NN = [68, 70, 72, 74, 76, 78, 80]
    nmaxs = [1000, 1000, 1000, 200, 200, 200, 200]
    for i, N in enumerate(NN):
        nmax = nmaxs[i]
        cmd = 'oaconference -N %d --ctype 2 --itype 3 --j1zero 1 --j3zero 1 -o dconferencej1j3-r --nmax %d' % (N, nmax)
        if verbose>=2:
            print(cmd)
        j = job(cmd, shorttag='N%d' % N)
        jj += [j]

    jobfile=join(scriptdir, 'subs.sh')
    with open(jobfile, 'wt') as fid:
        for i, j in enumerate(jj):
            outfile, s = createJobScript(j, index=i, verbose=0, queue='q72h', scriptdir=scriptdir)
            _ = fid.write('%s\n' % s)
    print('generated job file %s' % jobfile)


#%% Check

N = 26
kk = N - 1
cfile0 = 'cdesign-diagonal-%d-%d.oa' % (N, kk)
ee = oapackage.readarrayfile(os.path.join(outputdir, cfile0))
maxzpos = -1
ct = oapackage.conference_t(N, N, 0)
for ii, al in enumerate(ee):
    x = oapackage.extend_conference_matrix(al, ct, N - 1, 1, maxzpos)
    print('array %d: %d extensions' % (ii, x.nExtensions()))

#%% Generate example for C(8, 4)

if platform.node() == 'marmot':
    paperdir = '/home/eendebakpt/misc/oa/article-conference/'

    N = 8
    k = 4
    ll = oapackage.readarrayfile(
        os.path.join(outputdir, 'cdesign-%d-%d.oa' % (N, k)))

    l1 = oapackage.array2latex(
        np.array(ll[0]), header=1, comment='array 0 in C(8,4)', mode='pmatrix')
    l2 = oapackage.array2latex(np.array(ll[1]), header=1, mode='pmatrix')
    ss = r'\begin{figure}' + os.linesep
    ss += r'\begin{align*}' + os.linesep
    ss += 'X_0 = ' + l1  # + os.linesep
    ss += ', \quad' + os.linesep
    ss += 'X_1 = ' + l2  # + os.linesep
    ss += r'\end{align*}' + os.linesep
    ss += r'\caption{Designs in $C(8, 4)$}' + os.linesep
    ss += r'\label{figure:C84}' + os.linesep
    ss += r'\end{figure}' + os.linesep

    with open(os.path.join(paperdir, 'example-C84.tex'), 'w') as fid:
        _ = fid.write(ss)

    mz1 = oapackage.maxz(ll[0])
    print('array 1: maxz %d' % mz1)
    


