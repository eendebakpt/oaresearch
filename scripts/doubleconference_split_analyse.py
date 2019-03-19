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
import copy
import json_tricks as json

# setup data locations
resultsdir = join(os.path.expanduser('~'), 'oatmp')

# setup data locations
if 'VSC_SCRATCH' in os.environ.keys():
    vsccluster = 1
    print('we are running on the VSC cluster!')
    resultsdir = os.path.join(os.environ['VSC_DATA'], "oatmp")
else:
    vsccluster = 0

if not os.path.isdir(resultsdir):
    print('  resultsdir %s does not exist ... exiting' % resultsdir)
    sys.exit()


import oapackage
from oapackage import oahelper, splitDir, splitFile, splitTag  # reload(oahelper)
import oaresearch.research_conference

r = oapackage.log_print(-oapackage.SYSTEM, '')

# %% Helper functions


# %%
dobigcase = 48  # by default run a small case to test the scripts

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="output level", default=1)
parser.add_argument("-b", "--basedir", type=str, help="base directory")
parser.add_argument("--scriptdir", type=str,
                    default=None, help="base directory")
parser.add_argument("-s", "--statistics", type=int, default=0,
                    help="calculate statistics of generation")
parser.add_argument("-N", "--N", type=int, default=dobigcase,
                    help="number of splitting")
parser.add_argument("-ii", "--ii", type=int, default=6, help="only run part of top level")
parser.add_argument("-jj", "--jj", type=int, default=-1,
                    help="only run part of top+1 level")
parser.add_argument(
    "-js", "--jjstart", type=int, default=0, help="starting value of jj range")

parser.add_argument("--split", type=int, default=1,
                    help="use split mode in oa_depth_extend")
parser.add_argument("-X", "--X", type=int, default=0, help="debugging...")
parser.add_argument("--lastlevel", type=int, default=3, help="lastlevel")
parser.add_argument("--queue", type=str, default=None, help="default queue")
parser.add_argument("-c", "--ncores", type=int, default=0,
                    help="number of cores for job files")
parser.add_argument("--dozipsub", type=int, default=True, help="dozipsub")


#parser.add_argument('lvls', metavar='lvls', type=int, nargs='*',  help='list of split levels')
args = parser.parse_args()

if args.basedir is not None:
    resultsdir = args.basedir
    #print(args.basedir); exit(0)

N = args.N
print('double conference: case: %d, resultsdir %s' % (N, resultsdir))

# %%


def rfile(lvls, N, k, basetag='dconference'):
    if len(lvls) > 0:
        return basetag + '-' + splitFile(lvls[:]) + '-%d-%d' % (N, k) + '.oa'
    else:
        return basetag + '-%d-%d' % (N, k) + '.oa'


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


alljobs = []
cache = 1
verbose = 1
basetag = 'dconference'

outputdir = oapackage.mkdirc(os.path.join(
    resultsdir, 'doubleconference-%d' % (N, )))

print('--- Starting: case %d, outputdir %s' % (N, outputdir))


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


def gather_results(lvls, splitdata, paretofunction, verbose=1):
    level = len(lvls)
    k = splitdata[level]['n']
    if verbose:
        print('gather_results: level %s: getting %d subresults' % (lvls, k))
    splitdata[level]


if 0:
    print('TODO: gather results: maybe in C for efficiency??')
    gather_results([1], splitdata, None, verbose=1)

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


ll = listFiles(splitdata, args.ii, verbose=1)
os.chdir(outputdir)
oapackage.oainfo(os.path.join(ll[0]))

if 0:
    ll = listFiles(splitdata, 7, verbose=1)
    oapackage.oainfo(os.path.join(outputdir, ll[0]))
    arrays = oapackage.readarrayfile(os.path.join(outputdir, ll[0]))
    ll = listFiles(splitdata, 8, verbose=1)
    oapackage.oainfo(os.path.join(outputdir, ll[0]))

    ll = listFiles(splitdata, 12, verbose=1)
    oapackage.oainfo(os.path.join(outputdir, ll[0]))

# %% Generic Pareto/counting scheme


def calc_stats(ll, func, outputdir, verbose=1):
    """ Calculate statistics over generated designs

    Args:
        ll (list): list of files with designs
        func (function): function to apply to a list of designs
    """
    rr = []
    for ii, l in enumerate(ll):
        rverbose = 0
        lst = oapackage.readarrayfile(os.path.join(outputdir, l), rverbose)
        rr.append(func(lst))
        if verbose >= 2 or (verbose and ii % 20 == 0):
            print('count_stats: file %d/%d' % (ii, len(ll)))
            print('count_stats: %s' % (rr[-1]))
    return rr


import oaresearch


reload(oaresearch.research_conference)

from oaresearch.research_conference import SingleConferenceParetoCombiner


cache_dir = oapackage.mkdirc(os.path.join(outputdir, 'sc_pareto_cache'))

pareto_calculator = SingleConferenceParetoCombiner(outputdir, cache_dir=cache_dir, cache=True)
self = pareto_calculator


def evenodd_count(lst):
    """ Return number of foldover and number of even-odd designs """
    v = np.array([oapackage.isConferenceFoldover(al) for al in lst])
    return np.array([np.sum(v == True), np.sum(v == False)])


if 1:
    # ii=0
    #lst=oapackage.readarrayfile(os.path.join(outputdir, ll[ii]))

    verbose = 2
    rr = {}
    for k in range(2, splitdata['kmax'] + 1):
        ll = listFiles(splitdata, k, verbose=verbose)

        if 0:
            r = calc_stats(ll, evenodd_count, outputdir, verbose=verbose >= 2)
            totals = np.sum(np.array(r), axis=0)
            rr[k] = totals.tolist()
            print(Fore.BLUE + 'double conference N %d k %d: foldover, evenodd %s' % (N, k, totals,) + Fore.RESET)

        pareto_calculator.pre_calculate(ll)
        results = pareto_calculator.calculate(ll)
        pareto_calculator.write_combined_results(k, results)

        if results['N'] is not None:
            print(Fore.BLUE + 'k %d: Pareto arrays %d, Pareto classes %d' %
                  (k, results['npareto'], results['nclasses']) + Fore.RESET)


# %%

k = 7

ll = listFiles(splitdata, k, verbose=verbose)

r = calc_stats(ll, evenodd_count, outputdir, verbose=verbose >= 2)
