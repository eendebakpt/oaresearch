# -*- coding: utf-8 -*-
"""

Example script for calculating double conference matrices with split jobs

Pieter Eendebak <pieter.eendebak@gmail.com>
"""

# %% Load necessary packages
import os
import platform
import sys
import numpy as np
import time
from imp import reload
from os.path import join
import pdb
import argparse
import tempfile
from colorama import Fore
import itertools

# setup data locations
oadir = os.path.join(os.path.expanduser('~'), 'misc/oa/oacode/')
if platform.system() == 'Windows':
    oadir = os.path.join(os.path.expanduser('~'), 'svn', 'oapackage')

resultsdir = join(os.path.expanduser('~'), 'oatmp')

# setup data locations
if 'VSC_SCRATCH' in os.environ.keys():
    vsccluster = 1
    print('we are running on the VSC cluster!')
    oadir = os.path.join(os.environ['VSC_SCRATCH'], "OA-develop-source")
    resultsdir = os.path.join(os.environ['VSC_DATA'], "oatmp")
else:
    vsccluster = 0

if not os.path.isdir(oadir):
    print('  oadir %s does not exist ... exiting' % oadir)
    sys.exit()
sys.path.append(os.path.join(oadir, 'pythondevelop'))

if not os.path.isdir(resultsdir):
    print('  resultsdir %s does not exist ... exiting' % resultsdir)
    sys.exit()


import oapackage
print('oapackage: %s: %s' % (oapackage, oapackage.version()))

from oapackage import oahelper, splitDir, splitFile  # reload(oahelper)
from oapackage import splitTag, splitDir, splitFile
import researchOA
from researchOA import splitBase, job, numbersFile, gatherFilesList, gatherResults, checkLevel, gzipOA, doSplitFile
from researchOA import jobStatus, createJobScript
from researchOA import makeJobList


r = oapackage.log_print(-oapackage.SYSTEM, '')

# %% Helper functions


# %%
dobigcase = 20  # by default run a small case to test the scripts

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="output level", default=1)
parser.add_argument("-b", "--basedir", type=str, help="base directory")
parser.add_argument("--scriptdir", type=str,
                    default=None, help="base directory")
parser.add_argument("-s", "--statistics", type=int, default=0,
                    help="calculate statistics of generation")
parser.add_argument("-N", "--N", type=int, default=dobigcase,
                    help="number of splitting")
parser.add_argument("-ii", "--ii", type=int, default=- 1, help="only run part of top level")
parser.add_argument("-jj", "--jj", type=int, default=-1,
                    help="only run part of top+1 level")
parser.add_argument(
    "-js", "--jjstart", type=int, default=0, help="starting value of jj range")

parser.add_argument("--split", type=int, default=1,
                    help="use split mode in oa_depth_extend")
parser.add_argument("--lastlevel", type=int, default=3, help="lastlevel")
parser.add_argument("--queue", type=str, default=None, help="default queue")
parser.add_argument("-c", "--ncores", type=int, default=0,
                    help="number of cores for job files")


#parser.add_argument('lvls', metavar='lvls', type=int, nargs='*',  help='list of split levels')
args = parser.parse_args()

if args.basedir is not None:
    resultsdir = args.basedir
    #print(args.basedir); exit(0)

N = args.N
print('double conference: case: %d, resultsdir %s' % (N, resultsdir))

# %%


def compressOA(directory):
    lst = oapackage.findfiles(directory, '.*oa')
    for afile in lst:
        oapackage.oahelper.compressOAfile(os.path.join(directory, afile), verbose=0)


def rfile(lvls, N, k, basetag='dconference'):
    if len(lvls) > 0:
        return basetag + '-' + splitFile(lvls[:]) + '-%d-%d' % (N, k) + '.oa'
    else:
        return basetag + '-%d-%d' % (N, k) + '.oa'


def make_extend(lvls, splitdata, dosplit=False, verbose=1, makelog=True):
    N = splitdata['N']
    tag = splitTag(lvls)
    level = len(lvls)
    if len(lvls) == 0:
        k = splitdata[0]['k']
        km = splitdata[level]['k']
        startfile0 = 'dconference-%d-%d.oa' % (N, k)
    else:
        if level in splitdata:
            k = splitdata[level]['k']
        else:
            k = splitdata['kmax']
        km = k

        kprev = splitdata[level - 1]['k']
        startfile0 = 'dconference-' + splitFile(lvls) + '.oa'

    splitdir = splitDir(lvls[:-1])
    splitdirx = splitDir(lvls)
    startfile = os.path.join(outputdir, splitdir, startfile0)
    # splittag=splitTag(lvls)
    checkfile0 = rfile(lvls, N, k, basetag=basetag)
    checkfile = os.path.join(outputdir, splitdirx, checkfile0)

    checkfileZ0 = rfile(lvls, N, kprev + 1, basetag=basetag)
    checkfileZ = os.path.join(outputdir, splitdirx, checkfileZ0)

    oapackage.mkdirc(os.path.join(outputdir, splitdir))
    #dosplit = len(lvls) in splitdata
    spl = splitFile(lvls)
    basecmd = 'cd %s; ' % outputdir + 'mkdir -p %s;' % os.path.join(outputdir, splitdirx) + 'oaconference -N %d -k %d -f AB --ctype 2 --itype 3  --j1zero 1 --j3zero 1 -i %s -o %s --select 3 ' % (
        N, km, startfile, os.path.join(splitdirx, 'dconference-' + spl))

    if makelog:
        logfile = os.path.join(scriptdir, 'log-extend-%s.txt' % tag)
        basecmd += ' | tee %s' % logfile
    if verbose:
        print('make_extend %s: startfile %s' % (tag, startfile))
    if verbose >= 2:
        print('make_extend: cmd %s' % basecmd)
    cmd = basecmd
    # oapackage.runcommand(basecmd)
    if verbose:
        print('make_extend %s: checkfile %s' % (tag, checkfile0))
    if dosplit:
        print('   splitting resulting file')
        splitcmd, splitfile, checkfile = split_file([0], splitdata)
        # print(splitfile)
        # oapackage.runcommand(splitcmd)
        cmd += '; ' + splitcmd
        checkfile = os.path.join(outputdir, splitdirx, splitFile(lvls) + '.oa')
    return cmd, checkfile, startfile, {'basecmd': basecmd, 'checkfileZ': checkfileZ}


if 0:
    lvls = [19]
    # lvls=[]
    cmd, checkfile, startfile, r = make_extend(lvls, splitdata)


def split_file(lvls, splitdata, verbose=1):
    """ Split file at specified level """
    level = len(lvls)
    splitdir = splitDir(lvls)
    N = splitdata['N']
    if level in splitdata:
        k = splitdata[level]['k']
    else:
        k = splitdata['kmax']
    splitfile0 = os.path.join(
        splitDir(lvls), basetag + '-' + splitFile(lvls) + '-%d-%d' % (N, k) + '.oa')
    splitfile = os.path.join(outputdir, splitfile0)

    nsplit = splitdata[level]['n']
    spl = splitFile(lvls)
    sdir = os.path.join(outputdir, splitdir)
    splitcmd = 'cd %s; ' % outputdir + 'oasplit -n %d -f B -i %s -o %s' % (
        nsplit, splitfile, os.path.join(sdir, basetag + '-' + spl + '-sp%d' % level))
    if verbose:
        print('split_file: split %s into %d files' % (splitfile0, nsplit))
    checkfile = os.path.join(outputdir, sdir, basetag
                             + '-' + splitFile(lvls + [nsplit - 1]) + '.oa')
    return splitcmd, checkfile, splitfile, {'splitfile': splitfile}


if 0:
    splitcmd, checkfile, startfile, r = split_file([0], splitdata)
    print(splitcmd)
    print(checkfile)
    # oapackage.runcommand(splitcmd)

# %%
cache = 1
verbose = args.verbose
time.sleep(.1)

if N == 36:
    splitdata = {0: {'k': 5, 'n': 10}, 1: {
        'k': 7, 'n': 4}, 'N': N, 'kmax': N // 2}
elif N == 32 or N == 24 or N == 16 or N == 20:
    splitdata = {0: {'k': 4, 'n': 8}, 1: {
        'k': 7, 'n': 4}, 'N': N, 'kmax': N // 2}
elif N == 48:
    # first split at k=8 is about 4 hours of calculation, k=7 is pretty fast
    splitdata = {0: {'k': 7, 'n': 40}, 1: {
        'k': 9, 'n': 20}, 'N': N, 'kmax': N // 2}
else:
    raise Exception('splitdata not yet defined for N %d' % N)


# %%
homedir = os.getenv('HOME')
if args.scriptdir is None:
    if vsccluster:
        scriptdir = oapackage.mkdirc(
            join(resultsdir, 'doubleconference-%d-scripts' % N))
    else:
        scriptdir = oapackage.mkdirc(
            join(tempfile.mkdtemp(prefix='doubleconference-%d-' % N)))
else:
    scriptdir = oapackage.mkdirc(join(homedir, args.scriptdir))

print('--- Starting: case %d, ' %
      (N, ))

alljobs = []
cache = 1
verbose = 1
basetag = 'dconference'


# %% """ Step 1: extend to 5 columns, Even-odd """

print('--- Step 1 (sequential) ---')
sys.stdout.flush()
# create root and extend
kstart = 5
k = splitdata['kmax']
ctype = oapackage.conference_t(N, splitdata['kmax'], 1)
outputdir = oapackage.mkdirc(os.path.join(
    resultsdir, 'doubleconference-%d' % (N, )))
allcmdlogfile = os.path.join(outputdir, 'cluster-command-log.txt')
job.logfile = allcmdlogfile

# %% Initial stage

km = splitdata[0]['k']
basecmd = 'cd %s; ' % outputdir + \
    'oaconference -N %d -k %d --ctype 2 --itype 3 -f A --j1zero 1 --j3zero 1 -o %s --select 3 ' % (
        N, km, 'dconference')
checkfile = os.path.join(outputdir, 'dconference-%d-%d.oa' % (N, km))
splitcmd = 'cd %s; ' % outputdir + \
    'oasplit -n %d -f D -i %s -o %s' % (
        splitdata[0]['n'], checkfile, 'dconference-sp0')

if verbose >= 2:
    print(basecmd)
    print(splitcmd)


# %%
print(Fore.BLUE + 'compute level 1' + Fore.RESET)

if not oapackage.file_exists(checkfile) and cache:
    os.system(basecmd)
    print(Fore.BLUE + splitcmd + Fore.RESET)
    os.system(splitcmd)
    os.listdir(outputdir)
    print(Fore.BLUE + 'initial extension is computed' + Fore.RESET)
else:
    print(Fore.BLUE + 'initial extension is complete' + Fore.RESET)

#checkfile = ...

if not oapackage.file_exists(checkfile):
    raise Exception('base extention did not work')

if 0:
    lvls1 = [1, 208]
    joblistg1 = gatherResults(
        lvls1, outputdir, splitdata, adata=adata, verbose=2)
    j = joblistg1[0]
    print(j)
    j.analyse(verbose=2)

    exit(0)


# %%
print(Fore.BLUE + 'compute level 1' + Fore.RESET)

for ii in range(splitdata[0]['n']):
    lvls = [ii]
    tag = splitTag(lvls)
    cmd, checkfile, startfile, r = make_extend(lvls, splitdata, verbose=0)
    if oapackage.checkFilesOA(startfile) and not oapackage.checkFilesOA(checkfile):
        print('make job for extend of %s' % (lvls, ))
    # create job file
    j = job(cmd, jobtype='extend %s' % tag, checkfiles=[
            checkfile], checkfilesstart=[startfile])
    alljobs += [j]

    splitcmd, checkfile, startfile, r = split_file([ii], splitdata, verbose=0)
    j = job(splitcmd, jobtype='split %s' %
            tag, checkfiles=[checkfile], checkfilesstart=[startfile])
    if j.canrun() and not j.complete():
        print('make job for split of %s' % (lvls, ))
    alljobs += [j]

    cmd = researchOA.gzipOA(os.path.join(outputdir, splitDir(lvls)), unzip_text=True, cmdverbose=True)
    jc = job(cmd, jobtype='compress %s' % tag, checkfiles=[
        checkfile + '.gz'], checkfilesstart=[checkfile])
    jc.compressjob = True

    alljobs += [jc]

# %%
print(Fore.BLUE + 'compute level 2' + Fore.RESET)

for ii in range(splitdata[0]['n']):
    if ii > 3:
        # continue
        pass

    print('checking jobs for [%d, *]' % ii)
    for jj in range(splitdata[1]['n']):
        lvls = [ii, jj]
        level = len(lvls)
        tag = splitTag(lvls)
        cmd, checkfile, startfile, r = make_extend(lvls, splitdata, verbose=0)
        checkfileZ = r['checkfileZ']
        if oapackage.checkFilesOA(startfile) and not oapackage.checkFilesOA(checkfile):
            if verbose >= 2:
                print('make job for extend of %s' % (lvls, ))
        # create job file
        j = job(cmd, jobtype='extend %s' % tag, checkfiles=[
                checkfile], checkfilesstart=[startfile])
        alljobs += [j]

        cmd = researchOA.gzipOA(os.path.join(outputdir, splitDir(lvls)), unzip_text=True, cmdverbose=2)
        jc = job(cmd, jobtype='compress %s' % tag, checkfiles=[
            checkfileZ + '.gz'], checkfilesstart=[checkfile])
        jc.compressjob = True

        alljobs += [jc]
        # compressOA(os.path.join())

        if level in splitdata:
            splitcmd, checkfile, startfile, r = split_file(lvls, splitdata)
            j = job(splitcmd, jobtype='split %s' %
                    tag, checkfiles=[checkfile], checkfilesstart=[startfile])
            alljobs += [j]


# %%
if 0 and not vsccluster:
    for j in alljobs:
        j.execute = True
        time.sleep(0.05)
        j.runjob(cache=True)


# %%

def paretofunction(al):
    j4 = conferenceInvariant(al)
    f4 = al.FvaluesConference(jj=4)

    N = al.n_rows
    b4 = np.sum(np.array(j4)**2) / N**2
    #X2 = oapackage.array2secondorder(al)
    X2 = al.getModelMatrix(2)[:, (1 + al.n_columns):]

    r = np.linalg.matrix_rank(X2)

    if verbose >= 2:
        print('design %d: rank %d, b4 %.3f, F4 %s' % (ii, r, b4, f4))
    Q = np.array(al) * np.array(al)
    rx2q = np.linalg.matrix_rank(np.hstack((X2, Q)))

    presults.f4s += [f4]
    presults.ranks += [r]
    presults.b4s += [b4]

    presults.ranksX2Q += [rx2q]
    return presults


def gather_results(lvls, splitdata, paretofunction, verbose=1):
    level = len(lvls)
    k = splitdata[level]['n']
    if verbose:
        print('gather_results: level %s: getting %d subresults' % (lvls, k))
    splitdata[level]


# %%

def listFiles(splitdata, k, verbose=1):
    """ List all results files for a specified number of columns """
    for ii in range(10):
        if not ii in splitdata:
            level = ii
            break
        if k <= splitdata[ii]['k']:
            level = ii
            break
    nn = [splitdata[ii]['n'] for ii in range(level)]

    cc = list(itertools.product(*[range(i) for i in nn]))

    ll = [os.path.join(splitDir(c), rfile(c, N, k, basetag=basetag)) for c in cc]
    print('listFiles: k %d, level is %d, splits %s: %d file(s)' % (k, level, nn, len(ll)))
    return ll

#ll=listFiles(splitdata, 4, verbose=1)
#oapackage.oainfo(os.path.join(outputdir, ll[0]) )


if 0:
    ll = listFiles(splitdata, 7, verbose=1)
    oapackage.oainfo(os.path.join(outputdir, ll[0]))
    ll = listFiles(splitdata, 8, verbose=1)
    oapackage.oainfo(os.path.join(outputdir, ll[0]))

# %% Generic Pareto/counting scheme


def calc_stats(ll, func, verbose=1):
    rr = []
    for ii, l in enumerate(ll):
        if verbose:
            print('count_stats: %d/%d' % (ii, len(ll)))
        lst = oapackage.readarrayfile(os.path.join(outputdir, l))
        rr.append(func(lst))
        if verbose:
            print('count_stats: %s' % (rr[-1]))
    return rr


def evenodd_count(lst):
    v = np.array([oapackage.isConferenceFoldover(al) for al in lst])
    return np.array([np.sum(v == True), np.sum(v == False)])


if 0:
    # ii=0
    #lst=oapackage.readarrayfile(os.path.join(outputdir, ll[ii]))

    ll = listFiles(splitdata, 5, verbose=1)

    r = calc_stats(ll, evenodd_count)
    totals = np.sum(np.array(r), axis=0)
    print(totals)

# %% Load packages


jobs = alljobs


if args.queue is not None:
    queue = args.queue
else:
    queue = 'q1h'
    queue = 'q24h'
    queue = 'q72h'
    queue = None
    # queue='q7d'


if 0:
    for j in alljobs:
        if getattr(j, 'compressjob', False):
            print('enable compression job %s' % j)
            j.complete = lambda *args, **kwargs: False
            j.checkfiles = ['123abc_123']


jfile = makeJobList(scriptdir, jobs, ncores=args.ncores, queue=queue)

#os.system('. %s' % jfile)

pp = [j for j in alljobs if not j.complete() and j.canrun()]
#job.execute=True; [j.runjob() for j in alljobs]
#[j.runjob() for j in jobs]
