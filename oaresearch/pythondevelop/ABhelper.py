"""
Created on Mon Jun 13 14:12:51 2011

@author: eendebakpt
"""

""" Load necessary packages """


#x = np.seterr(invalid='raise')
# np.geterr()




import fileinput
import os
import sys
from os.path import basename
import matplotlib.pyplot as plt
import numpy as np
import oalib
import oapackage.oahelper as oahelper
from oapackage.oahelper import floatformat, checkFiles, getArrayFile
from oapackage.oahelper import *
def extendFile(afile, ostr, configfile, verbose=1, cache=1, logfile=None, cmdlog=None):
    """ Extend arrayfile with dynamic filtering """
    print('incomplete')
    if logfile == None:
        cmd = 'oaextendsingle -c %s -l 2 -r %s -o %s' % (
            configfile, afile, ostr)
    else:
        cmd = 'oaextendsingle -c %s -l 2 -r %s -o %s | tee %s' % (
            configfile, afile, ostr, logfile)
    nextfile = ''
    if verbose >= 2:
        print(cmd)
    if cmdlog != None:
        cmdlog.write('# Extend \n' + scriptCheck(cmd, nextfile) + '\n\n')
        cmdlog.flush()
    sys.stdout.flush()  # needed for tee command
    if not checkFiles(nextfile, cache=cache):
        os.system(cmd)
    else:
        if verbose >= 2:
            icmd = 'oainfo %s' % nextfile
            os.system(icmd)
    return cmd


def gma2str(gmadata, t=None, sformat=None):
    """ Convert GWLP value to string format """
    print('legacy function')
    if gmadata is None:
        return '-'
    gmadata[gmadata < 0] = 0
    if not(np.abs(gmadata[0] - 1) < 1e-12 and np.abs(gmadata[1]) < 1e-12):
        print('warning: data are not good GWPL data!!!!')
    bgma = np.around(gmadata, decimals=12)
    if not t is None:
        bgma = bgma[(t + 1):]
    if sformat is None:
        gstr = ','.join([floatformat(v, mind=2, maxd=4) for v in bgma])
    else:
        gstr = ','.join([sformat % v for v in bgma])
    return gstr


def analyseArraysCache(afile, verbose=1):
    anafile = analyseFile(afile)

    data = loadAnalysisFile(anafile)
    a = data[:, 1]
    return a


def calculateAthreshold(Afinal, kfinal, k, L=1):
    """ Calculate D-eff threshold at specified number of columns """
    mfinal = 1 + kfinal + kfinal * (kfinal - 1) / 2
    m = 1 + k + k * (k - 1) / 2
    Cfinal = Afinal**mfinal
    Cthr = Cfinal / (L**(kfinal - k))
    Athr = Cthr**(float(1) / float(m))
    return Athr


def drawAvalues(aval, fig=1):
    """ From a values in plot """
    fig = plt.figure(fig)
    plt.clf()
    if aval.size > 0:
        plt.plot(aval, '.b', label='a')
    else:
        plt.plot([], '.b', label='a')
    plt.ylabel('D-efficiency', fontsize=15)
    plt.xlabel('Index', fontsize=15)
    if aval.size > 0:
        plt.title('max D-eff value %.3f' % aval.max())
    else:
        plt.title('max D-eff value %.3f' % 0)
    plt.draw()
    return fig


def Dvalue2col(A, k, kout, Lmax=1):
    """ Convert D-efficiency to equivalent value at another number of columns """
    m = 1 + k + k * (k - 1) / 2
    mout = 1 + kout + kout * (kout - 1) / 2
    dc = kout - k
    Aout = (A**(float(m)) * (Lmax**dc))**(1. / mout)

    return Aout


def nextDthr(A, k):
    m = 1 + k + k * (k - 1) / 2
    kn = k + 1
    mn = 1 + kn + kn * (kn - 1) / 2
    An = A**(float(m) / mn)
    return An


def plotAthresholdsY(Afinal, kfinal, k=None, yl=None, label='Threshold for final D-efficiency value'):
    xl = plt.xlim()
    yl = plt.ylim()
    if k == None:
        k = kfinal
    Afinaln = calculateAthreshold(Afinal, kfinal, k)
    ph = plt.plot(xl, [Afinaln, Afinaln], '--r', linewidth=3, label=label)
    return ph


def plotAthresholds(Afinal, kfinal, k=None, yl=None):
    """ Plot thresholds for AB figure """
    if yl == None:
        yl = plt.ylim()
    plt.plot([Afinal, Afinal], yl, '-g')
    Acurr = Afinal
    if not k == None:
        mfinal = 1 + kfinal + kfinal * (kfinal - 1) / 2
        m = 1 + k + k * (k - 1) / 2
        Acurr = Afinal**(float(mfinal) / m)
        plt.plot([Acurr, Acurr], yl, '--r', label='Threshold for D-efficiency')
    return Acurr


def parseProcessingTimeOld(logfile, verbose=0):
    if verbose:
        print('ERROR: do not use this functioN!')
    fileinput.close()
    tstart = None
    tend = None
    for line in fileinput.input([logfile]):
        if line.startswith(' TIME'):
            if verbose:
                print(line)
            tstart = float(line.split(' ')[2])
            break
    dtt = tstart
    return dtt


def setABfigure(figid=None, fontsize=13):
    """ Set labels for AB figure """
    if not figid == None:
        plt.figure(figid)
    plt.xlabel('D-efficiency', fontsize=fontsize)
    plt.ylabel('log(average variation inflation factor)', fontsize=fontsize)
    plt.title('Scatterplot', fontsize=16)
    plt.xlim([0, 1.1])


def plotABfigure(a, b, figid=1, verbose=1, fontsize=13):
    """ Plot """
    plt.figure(figid)
    plt.clf()
    gidx = b > 0
    plt.plot(a[gidx], np.log(b[gidx]), '.')
    setABfigure(figid, fontsize=fontsize)


def plotAsAboundaries(k):
    alpha = (1 + (k - 1) + (k - 1) * (k - 2) / 2) / (1 + k + k * (k - 1) / 2)
    xx = np.arange(0, 1, .01)
    yy = xx**alpha
    h = plt.plot(xx, yy, '--r', label='Boundary')
    return h


def plotABboundaries(Acc=None):
    """ Plot boundaries for D-eff and A-eff values """
#    xl=plt.xlim();
    xl = [0, 1]
    yl = plt.ylim()
    plt.plot([1, 1], yl, '--r')
    dstep = .005
    xx = np.arange(dstep, 1, dstep)
    yy = 1 / xx
    idx = np.log(yy) < yl[1]
    hborder = plt.plot(xx[idx], np.log(yy[idx]), '--r', label='Boundary')

    if not Acc == None:
        hacc = plt.plot(
            [Acc, Acc], yl, '--m', label='Numerical accuracy limit')


def ABcalclist(sols, verbose=0):
    nn = sols.size()
    a = np.zeros(nn)
    b = np.zeros(nn)
    for ii, X in enumerate(sols):
        if verbose:
            print('ABcalclist: %d' % ii)
        (A, B) = ABcalc(X.getarray())
        a[ii] = A
        b[ii] = B
    return a, b


def ABcalc(X, verbose=0):
    """ Calculate D-efficiency and A-efficiency values of an array """
    Y = X.copy()
    n = Y.shape[0]
    k = Y.shape[1]
    Y = Y.astype(float)
    Y -= .5
    Y *= 2
    ntf = k * (k - 1) / 2
    tf = np.zeros((n, ntf))
    kk = 0
    for ii in range(0, k):
        for jj in range(0, ii):
            tf[:, kk] = Y[:, ii] * Y[:, jj]
            kk = kk + 1
    c = np.ones((n, 1))
    x = np.hstack((c, Y, tf))
    m = x.shape[1]

    (u, s, v) = np.linalg.svd(x)
    if np.any(s < 1e-15):
        # singular matrix
        A = 0
        B = 0
    else:
        A = np.exp(np.sum(2 * np.log(s)) / m - np.log(n))
        B = n * np.sum((1 / s)**2) / m

    if 0:
        xtx = np.asmatrix(x.transpose()) * np.asmatrix(x)
        # old method
        if verbose:
            print('xtx')
            print(xtx)
        tmp = np.linalg.det(xtx / 32)

        if verbose >= 2:
            print('tmp %e' % tmp)
        if np.abs(tmp) < 1e-20:
            if verbose:
                print('tmp %f' % tmp)
            tmp = 0
        A = tmp**(1. / m)
        if verbose >= 2:
            print('A %em m %d' % (A, m))
        if tmp != 0:
            B = n * xtx.I.diagonal().sum()
        else:
            B = 0
    if verbose:
        print('D-eff %f A-eff %f' % (A, B))
    return (A, B)


def x2xf(X, convert=0):
    """ Apply intercept and second order interactions to design """
    Y = X.copy()
    n = Y.shape[0]
    k = Y.shape[1]
    if convert:
        Y = Y.astype(float)
        Y -= .5
        Y *= 2
    ntf = k * (k - 1) / 2
    tf = np.zeros((n, ntf))
    kk = 0
    for ii in range(0, k):
        for jj in range(0, ii):
            tf[:, kk] = Y[:, ii] * Y[:, jj]
            kk = kk + 1
    c = np.ones((n, 1))
    x = np.hstack((c, Y, tf))
    return x


def processBlocks(nblocks, outfiles, analysisfiles, verbose, cverbose, cache=1, cmdfile=None):
    """ Process data in blocks """

    if cmdfile is not None:
        cfid = open(cmdfile, 'wt')

    for nn in range(0, nblocks):
        if verbose:
            print('processBlocks: block %d' % nn)
        afile = 'base-split-%d' % nn + '.oa'
        outfile = outfiles[nn]
        anafile = analysisfiles[nn]

        outfiler = getArrayFile(outfile)

        cmd = 'echo "block %d"; nice oaextendsingle -f B -r %s -o res%d -l 2 | tee log-base%d.txt' % (
            nn, afile, nn, nn)
        if cache and (os.path.exists(outfile) or os.path.exists(outfile + '.gz')):
            if verbose >= 2:
                print('   output %s exists' % outfiler)
            if cverbose >= 2:
                print(cmd)
        else:
            if cmdfile is not None:
                cfid.write(cmd + '\n')
            if nn % 2 != 0:
                cmd = cmd + ' &'
            if cverbose >= 2:
                print(cmd)

            continue
        cmd = 'echo "analyse block %d"; nice oaanalyse -j -1 -r -a ana-%d %s | tee log-base-analysis-%d.txt' % (
            nn, nn, outfiler, nn)
        if os.path.exists(anafile) and cache and nn > 0:
            if verbose >= 2:
                print('   analysis file %s exists' % anafile)
            if cverbose >= 2:
                print(cmd)
        else:
            if cmdfile is not None:
                cfid.write(cmd + '\n')
            if nn % 2 != 0:
                cmd = cmd + ' &'
            if cverbose:
                print(cmd)
            continue


def splitFunc(X, verbose=0):
    """ Variation of x2xf """
    Y = X.copy()
    n = Y.shape[0]
    k = Y.shape[1]
    Y = Y.astype(float)
    #m = 1 + k + k * (k - 1) / 2
    Y -= .5
    Y *= 2
    c = np.ones((n, 1))
    ntf = k * (k - 1) / 2
    X2 = np.zeros((n, ntf - (k - 1)))
    kk = 0
    for ii in range(0, k - 1):
        for jj in range(0, ii):
            X2[:, kk] = Y[:, ii] * Y[:, jj]
            kk = kk + 1
    E = Y[:, (k - 1):k]
    E2 = np.zeros((n, (k - 1)))
    kk = 0
    for ii in [k - 1]:
        for jj in range(0, ii):
            E2[:, kk] = Y[:, ii] * Y[:, jj]
            kk = kk + 1
    x = np.hstack((c, Y[:, 0:-1], X2, E, E2))
    return x


def parseABarrayfull(X, subcols=None, verbose=1):
    """ Parse an array, return a and asub values """
    A, B = ABcalc(X.copy())
    nc = X.shape[1]
    asub = np.zeros(nc)
    bsub = np.zeros(nc)
    if subcols == None:
        subcols = range(0, nc)
    for ii in subcols:
        Xs = np.delete(X, np.s_[ii], axis=1)
        As, Bs = ABcalc(Xs)
        asub[ii] = As
        bsub[ii] = Bs
    return(A, B, asub, bsub)


def parseABarray(af, nn, verbose=1):
    """ Parse a list of arrays (from arrayfile) """
    al = oalib.array_link(af.nrows, af.ncols, 0)
    a = np.zeros((nn, 1))
    b = np.zeros((nn, 1))
    asub = np.zeros((nn, 1))
    bsub = np.zeros((nn, 1))
    verbose = 1
    for ii in range(0, nn):
        if verbose >= 1:
            if ii % 400 == 0:
                print('fraction %.3f' % (float(ii) / float(nn)))
        t = af.read_array(al)
        X = al.getarray()
        A, B = ABcalc(X.copy())
        As, Bs = ABcalc(X[:, 0:-1])
        if A > 0:
            if verbose >= 2:
                print('array %d: D-efficiency %.3f A-efficiency %.3f' %
                      (ii, A, B))
        a[ii] = A
        b[ii] = B
        asub[ii] = As
        bsub[ii] = Bs
    return(a, b, asub, bsub)
