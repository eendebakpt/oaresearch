"""

Example script for calculating even-odd arrays for OA(N, 3, 2^a)

The script consists of 5 steps:

    1. Extend to OA(N; 3; 2^5), discard the arrays that have J5==0
    2. Extend the selected arrays to OA(N; 3; kinitial)
    2b. Split the resulting arrays in blocks
    3. For each block available:
        3a. Extend to OA(N; 3; 2^knext) using the algorithm ORDER_J5X
        3b. Split the resulting arrays into blocks
    4. For each block available:
        4a. Extend to OA(N; 3; 2^klevel) using the algorithm ORDER_J5X and a special algorithm (depth mode)
        *** OR ***
        4a.i Extend to OA(N; 3; 2^(kinitial+delta) )
        4a.ii Split into a number of blocks
        4a.iii Extend each of the blocks
        4a.iv Merge extended blocks into original format

    5. Collect the results

Pieter Eendebak <pieter.eendebak@gmail.com>
"""

# %% Load necessary packages

import argparse
import os
import platform
import sys
import time
from importlib import reload
from os.path import join

import numpy as np
import oapackage
from oapackage import oahelper  # reload(oahelper)
from oapackage import MODE_J5ORDERX

from oaresearch.pythondevelop.researchOA import makeJobList

import oaresearch.pythondevelop.researchOA as researchOA
from oaresearch.pythondevelop.functions_evenodd import generateLevelNext
from oaresearch.pythondevelop.researchOA import (
    analyseFile,
    checkFiles,
    checkLevel,
    combineNumbers,
    createJobScript,
    doSplitFile,
    evenoddAnalysePartialRun,
    evenoddAnalyseRun,
    makeJobList,
    evenoddCases,
    findfiles,
    gatherFilesList,
    gatherResults,
    gzipOA,
    job,
    jobStatus,
    loadAnalysisFile,
    mkdirc,
    numbersFile,
    parseParetoList,
    runcommand,
    splitBase,
    splitdir,
    splitFile,
    splitname,
    splitTag,
)

# %% setup data locations
resultsdir = join(os.path.expanduser("~"), "oatmp")

# setup data locations
if "VSC_SCRATCH" in os.environ.keys():
    vsccluster = 1
    print("we are running on the VSC cluster!")
    resultsdir = os.path.join(os.environ["VSC_SCRATCH"], "OA-run")
else:
    vsccluster = 0


if not os.path.isdir(resultsdir):
    print("  resultsdir %s does not exist ... exiting" % resultsdir)
    sys.exit()


r = oapackage.log_print(-oapackage.SYSTEM, "")

# %% Setup arguments

# dobigcase = 56  # by default run a small case to test the scripts
dobigcase = 48  # by default run a small case to test the scripts
dobigcase = 64

print("evenodd_cluster: command line usage evenodd_cluster.py -N [N]")

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="output level", default=1)
parser.add_argument("-b", "--basedir", type=str, help="base directory")
parser.add_argument("--scriptdir", type=str, default=None, help="base directory")
parser.add_argument("-s", "--statistics", type=int, default=0, help="calculate statistics of generation")
parser.add_argument("-N", "--N", type=int, default=dobigcase, help="number of splitting")
parser.add_argument("-ii", "--ii", type=int, default=-1, help="only run part of top level")
parser.add_argument("-jj", "--jj", type=int, default=-1, help="only run part of top+1 level")
parser.add_argument("-js", "--jjstart", type=int, default=0, help="starting value of jj range")

parser.add_argument(
    "--discardJ5", type=int, default=-1, help="discard arrays with J5=max after specified number of columns"
)
parser.add_argument("--split", type=int, default=1, help="use split mode in oa_depth_extend")
parser.add_argument("-X", "--X", type=int, default=0, help="debugging...")
parser.add_argument("--lastlevel", type=int, default=3, help="lastlevel")
parser.add_argument("--queue", type=str, default=None, help="default queue")
parser.add_argument("-c", "--ncores", type=int, default=0, help="number of cores for job files")
parser.add_argument("--dozipsub", type=int, default=True, help="dozipsub")

parser.add_argument("--execute", type=int, default=0, help="execute generated jobs")


# parser.add_argument('lvls', metavar='lvls', type=int, nargs='*',  help='list of split levels')
args = parser.parse_args()
args.execute = 1


if args.basedir is not None:
    resultsdir = args.basedir
    # print(args.basedir); exit(0)

dobigcase = args.N
print("evenodd: case: %d, resultsdir %s" % (dobigcase, resultsdir))

dryrun = 0
cache = 1
verbose = args.verbose
time.sleep(0.1)

# method used
base = "depth"

if 0:
    # hack
    compresssubdir = 1
    # job.execute=True

splitdata, iisel, jjsel = evenoddCases(dobigcase, strength=3, lastlevel=args.lastlevel)
N = splitdata["N"]
strength = splitdata["strength"]

if args.ii > -1:
    iisel = [args.ii]
if args.jj > -1:
    jjsel = list(range(args.jjstart, min(args.jj + 1, splitdata[1]["n"])))

paretomethod = oapackage.PARETOFUNCTION_J5
discardJ5 = args.discardJ5  # -1

# %%
if dobigcase == 32:
    cache = 0
    dryrun = 0
if dobigcase == 40:
    compresssubdir = 1
    cache = 1

if dobigcase == 64:

    if vsccluster:
        compresssubdir = 1
        dryrun = 1
    else:
        compresssubdir = 0
        dryrun = 0

if vsccluster:
    compresssubdir = 1
    dryrun0 = 1
    dryrun = 1
else:
    dryrun0 = 0
    compresssubdir = 1
    if dryrun == 1:
        compresssubdir = 1


# %%
homedir = os.getenv("HOME")
if args.scriptdir is None:
    scriptdir = oapackage.mkdirc(join(homedir, "b%d" % N))
else:
    scriptdir = oapackage.mkdirc(join(homedir, args.scriptdir))


# %% """ Even-odd """

print("--- Starting: case %d, dryrun %d, dryrun0 %d, compresssubdir %d" % (dobigcase, dryrun, dryrun0, compresssubdir))

alljobs = []


# %% """ Step 1: extend to 5 columns, Even-odd """

print("--- Step 1 (sequential) ---")
sys.stdout.flush()
# create root and extend
kstart = 5
k = splitdata["kmax"]
adata = oapackage.arraydata_t(oapackage.intVector([2] * k), N, strength, k)
oaoptions = oapackage.OAextend()
outputdir = oapackage.mkdirc(os.path.join(resultsdir, "eocluster-%s-%d-t%d" % (base, N, strength)))
allcmdlogfile = os.path.join(outputdir, "cluster-command-log.txt")
job.logfile = allcmdlogfile

with open(os.path.join(outputdir, "code.txt"), "wt") as fid:
    fid.write(f"#generated by evenodd_cluster.py\n")
    fid.write(f"# oapackage {oapackage.version()}\n")


# %%

if not os.path.exists(outputdir):
    os.mkdir(outputdir)
adata.writeConfigFile(os.path.join(outputdir, "oaconfig.txt"))
adatastart = oapackage.arraydata_t(adata, kstart)
adatastart.writeConfigFile(os.path.join(outputdir, "oaconfig-start.txt"))

# create file with root array
adata0 = oapackage.arraydata_t(adata, adata.strength)
al = adata.create_root()
afile = "result-%s" % adata0.idstr() + ".oa"
afile5 = "result-%s" % oapackage.arraydata_t(adata, 5).idstr() + ".oa"
if oapackage.checkFilesOA(os.path.join(outputdir, afile)):
    pass
else:
    oapackage.writearrayfile(os.path.join(outputdir, afile), al, oapackage.ABINARY)

cmdlogfile = "%s-extend.txt" % splitBase([])
cmd = "cd %s; oaextendsingle -f B -l 2 -m %d -c oaconfig-start.txt --maxk 5 | tee %s" % (
    outputdir,
    MODE_J5ORDERX,
    cmdlogfile,
)

j = job(cmd, jobtype="extend initial", checkfiles=[join(outputdir, afile5)], checkfilesstart=[join(outputdir, afile)])

os.chdir(outputdir)

if j.complete():
    print("  file %s already there" % afile)
else:
    j.execute = True
    time.sleep(0.05)
    j.runjob(cache=True)


afile = "result-%s" % adatastart.idstr() + ".oa"
eostartfile = "eo-%s" % oapackage.arraydata_t(adata, 5).idstr() + ".oa"

if oahelper.checkFiles(os.path.join(outputdir, eostartfile), cache=cache):
    if verbose:
        print("  eostartfile: %s already there" % eostartfile)
else:
    anafile = analyseFile(afile, method="gwlp", verbose=1, cache=cache)
    gwlp = loadAnalysisFile(anafile[0])

    ii = (np.abs(gwlp[:, -1]) > 1e-5).nonzero()[0]
    eo_indices = [int(x) for x in ii]
    #idxvec = oapackage.intVector(ww)
    sols0 = oapackage.arraylist_t()
    for ii in eo_indices:
        al = oapackage.oalib.selectArrays(afile, ii)
        sols0.push_back(al)
    oapackage.writearrayfile(eostartfile, sols0, oapackage.ABINARY)


# %%  Step 2. Extend the selected arrays to OA(N; 3; kinitial)

print('--- Step 2 (sequential) ---')
kinitial = splitdata['kinitial']
cmdlogfile = 'logstart.txt'

adatastart = oapackage.arraydata_t(adata, kinitial)
adatastart.writeConfigFile(os.path.join(outputdir, "oaconfig-start.txt"))
cmd = f"cd {outputdir}; oaextendsingle -r {eostartfile} -f D -l 2 -m {oapackage.MODE_J5ORDERX} -c oaconfig-start.txt -o eo | tee {cmdlogfile}"
afile = "eo-%s" % adatastart.idstr() + ".oa"
j = job(cmd=cmd, checkfiles=afile, checkfilesstart=eostartfile, jobtype="extend initial")
j.execute = True
res = j.runjob(cache=cache)

os.chdir(outputdir)
dt0 = oahelper.parseProcessingTime(cmdlogfile)
print("time initial extension: %.1f [s], %d arrays " % (dt0, oapackage.nArrays(afile)))

#%% Research:

startfile = "result-48.2-2-2-2-2.oa"
cmd = f"oaextendsingle -r {startfile} -f D -l 2 -c oaconfig.txt -o full"
print(cmd)

startfile = "result-48.2-2-2-2-2.oa"
cmd = f"oaextendsingle -r {startfile} -f D -l 2 -m {oapackage.MODE_J5ORDERX} -c oaconfig.txt -o full"
print(cmd)

# only even-odd
cmd = f"oaextendsingle -r {eostartfile} -f D -l 2 -m {oapackage.MODE_J5ORDERX} -c oaconfig.txt -o onlyeo"
print(cmd)

# %%
if 0:
    ll = oapackage.readarrayfile(eostartfile)
    ll = oapackage.readarrayfile(r"result-48.2-2-2-2-2.oa")  # [2, 5, 6, 8, 9]
    for idx, d in enumerate(ll):
        oapackage.writearrayfile(f"single{idx}.oa", ll[idx])
        cmd = f"oaextendsingle -r single{idx}.oa -f D -l 2 -m {oapackage.MODE_LMC_2LEVEL} -c oaconfig.txt -o single"
        print(cmd)
        cmd = f"oaextendsingle -r single{idx}.oa -f D -l 2 -m {oapackage.MODE_J5ORDERX}-c oaconfig.txt -o single-eo"
        print(cmd)

    # even-odd N=48
    # 0: 0.9 [s]
    # 1: 2.8 [s]
    # 2: 0.2 [s]
    # 3: 0.0 [s]
    # 4: 877.9 [s]

# %% """ Step 2b: split the initial extension """

print("--- Step 2b (sequential): split the initial extension ---")

level = 0

n = splitdata[0]["n"]
rfile = splitname([n - 1])
cmd = "oasplit  -v 1 -f D --nb 1 -i %s -n %d -o sp0" % (afile, n)
j = job(cmd=cmd, checkfiles=rfile, checkfilesstart=afile, jobtype="split level base")
j.execute = True
res = j.runjob(cache=cache)


# %% Step 3: extend with normal extension, then split """


print("--- Step 3 ---")

klevel = splitdata["klevel"]
level = 1
configfile = "oaconfig%d.txt" % klevel
adatax = oapackage.arraydata_t(adata, klevel)
adatax.writeConfigFile(os.path.join(outputdir, configfile))

alljobs = []
for ii in iisel:
    rfile = splitname([ii])
    if verbose:
        print("split level %d: file %s (extend and split)" % (level, rfile))
    ebase = rfile.replace(".oa", "-extend")
    edir0 = splitdir([ii])
    edir = splitdir([ii])
    tag = splitTag([ii])
    _ = mkdirc(edir)

    splitfile = "splitted-%s.zip" % edir0
    if oahelper.checkFiles(os.path.join(outputdir, edir0, splitfile), cache=cache):
        print("  split file %s in old format ! " % splitfile)
        continue
    splitfile = "splitted-%s.zip" % splitTag([ii])
    if oahelper.checkFiles(os.path.join(outputdir, edir0, splitfile), cache=cache):
        print("  split directory %s already calculated " % edir0)
        continue

    afile = os.path.join(edir, "%s-%s" % (ebase, adatax.idstr() + ".oa"))
    cmdlogfile = os.path.join(edir, rfile.replace(".oa", "-extend.txt"))
    cmd = f"cd {outputdir}; mkdir -p {edir}; oaextendsingle -f D -l 2 -m {oapackage.MODE_J5ORDERX} -c {configfile} -r {rfile} -o {os.path.join(edir, ebase)} | tee {cmdlogfile};\n"

    if 1:
        j = job(
            cmd=cmd,
            jobtype="extend %s" % tag,
            shorttag="E%s" % tag,
            checkfiles=[os.path.join(outputdir, afile)],
            checkfilesstart=[rfile],
        )
        # j.execute=True

        # print(j); print(j.execute); exit(0)

        j.runjob()
        alljobs += [j]

        lvls = [ii]
        jjs, _, _ = generateLevelNext(
            adata,
            lvls,
            splitdata,
            outputdir,
            allcmdlogfile=cmdlogfile,
            discardJ5=discardJ5,
            paretomethod=paretomethod,
            ncores=4,
            priorextend=False,
            verbose=1,
        )
        for jx in jjs:
            jx.queue = "q24h"
        [jx.runjob(verbose >= 2) for jx in jjs]

    alljobs += jjs


gjobs, jobs = jobStatus(alljobs)

# [j.runjob() for j in jobs]

# %%
if args.statistics:
    print("--- Step 3b (sequential): statistics of first split ---")

    ntotal = 0
    nfiles = 0
    totaltime = 0
    for ii in range(0, splitdata[0]["n"]):
        lvls = [ii]
        rfile = splitname([ii])
        ebase = rfile.replace(".oa", "-extend")
        edir = splitdir([ii])
        afile = os.path.join(edir, "%s-%s" % (ebase, adatax.idstr() + ".oa"))
        cmdlogfile = os.path.join(edir, rfile.replace(".oa", "-extend.txt"))
        dt0 = oahelper.parseProcessingTime(cmdlogfile)

        if dt0 >= 0:
            totaltime += dt0
        v = oapackage.nArrays(afile)
        if v >= 0:
            print("file %d: %d arrays (%d columns)" % (ii, v, adata.ncols))
            ntotal += v
            nfiles += 1
    if nfiles > 0:
        sfac = 1 / (float(nfiles) / splitdata[0]["n"])
        print(
            "stage 1: %d/%d files extended, # arrays %d (estimate %e)"
            % (nfiles, splitdata[0]["n"], ntotal, ntotal * sfac)
        )
        print("  time %.1f [s], estimate %.1f [h] " % (totaltime, float(totaltime) * sfac / 3600.0))
        print("  time estimated: cluster %.1f [d]" % (float(totaltime) * sfac / (100 * 24.0 * 3600.0)))


# %% Step 4: extend with special code """


# klevel2=12
print("--- Step 4 ---")


# configfile='oaconfig%d.txt' % splitdata['klevel2']
adatax2 = oapackage.arraydata_t(adata, splitdata["klevel2"])
adatax2idstr = adatax2.idstr()
# adatax.writeConfigFile(os.path.join(outputdir, configfile ))

# print(iisel); print(jjsel); raise

level = 2
joblist = []
for ii in iisel:
    # hack
    if ii % 40 == 0:
        print("calculating %d/%d" % (ii, splitdata[0]["n"]))

    edir0 = splitdir([ii])
    splitfile = "splitted-%s.zip" % edir0
    if oahelper.checkFiles(os.path.join(outputdir, edir0, splitfile), cache=cache):
        print("  split file %s in old format ! " % splitfile)
        continue
    splitfile = "splitted-%s.zip" % splitTag([ii])
    if oahelper.checkFiles(os.path.join(outputdir, edir0, splitfile), cache=cache):
        print("  split directory %s already calculated " % edir0)
        continue

    runcomplete = 1
    for jj in jjsel:
        if splitdata["lastlevel"] < 3:
            continue
        lvls2 = [ii, jj]
        tag2 = splitTag(lvls2)
        edir = splitdir(lvls2)
        sn = splitname(lvls2)
        cmdlogfile = os.path.join(edir, sn.replace(".oa", "-generate.txt"))

        # job.execute=True
        joblist, runcomplete2, cmd2 = generateLevelNext(
            adata,
            lvls2,
            splitdata,
            outputdir,
            allcmdlogfile=cmdlogfile,
            splitmode=args.split,
            ncores=8,
            priorextend=False,
            discardJ5=discardJ5,
        )
        alljobs += joblist

        # print('lvls %s: %s' % (lvls, joblist))
        # print('evenodd_cluster: xxxxxx after generate:' ); gjobs, jobs = jobStatus(alljobs);

        if 0:
            if N <= 40:
                _ = [jx.runjob() for jx in joblist]
        # checkLevel(lvls2, splitdata, adata, outputdir, verbose=1)

        legacyfile = join(outputdir, edir, "legacy.txt")
        if checkFiles(legacyfile):
            if verbose >= 2:
                print("generateLevelNext: %s: legacy format (generation complete)" % splitTag(lvls2))
            continue
        lockfile = join(outputdir, edir, "lockfile-%s.txt" % tag2)
        if checkFiles(lockfile):
            if verbose >= 1:
                print("  !! generateLevelNext: %s: lockfile exists!" % splitTag(lvls2))
            continue
        joblistg = gatherResults(lvls2, outputdir, splitdata, adata=adata, ncores=10, verbose=2)
        alljobs += joblistg
        # print('xxxxxx after gather:' ); gjobs, jobs = jobStatus(alljobs);

        print("lvls %s: %s" % (lvls, joblistg))
        joblistg[0].canrun(verbose=1)
        # [j.runjob() for j in joblistg]

        # checkLevel(lvls2, splitdata, adata, outputdir, verbose=1)

    print("  completed generation of block %d (runcomplete %d), collecting results" % (ii, runcomplete))
    # calculate numbers, select pareto optimal designs, compact directories

    kmax = adata.ncols
    kmin = klevel + 1

    # job.execute=True
    print("evenodd_cluster: %s:" % researchOA.splitTag([ii]))
    gjobs, jobs = jobStatus(alljobs)  # exit(0)

    lvls1 = [ii]
    joblistg1 = gatherResults(lvls1, outputdir, splitdata, adata=adata, verbose=2, ncores=4, dozip=args.dozipsub)
    alljobs += joblistg1
    # [j.runjob() for j in joblistg1]
    # print(joblistg1)

if N <= 40:
    _ = [jx.runjob() for jx in alljobs]

if verbose >= 2:
    gjobs, jobs = jobStatus(alljobs)

# [ jx.analyse(verbose=1) for jx in alljobs]
# print(iisel); print(jjsel)

# %% Step 5 (sequential): collect results """
print("--- Step 5 ---")

# Calculate pareto optimal arrays for initial sequence
for kk in range(5, splitdata["kinitial"] + 1):
    if dryrun0:
        continue
    print("pareto: %d columns" % kk)
    sys.stdout.flush()
    adata0 = oapackage.arraydata_t(adata, kk)
    xfile = "eo-%s.oa" % (adata0.idstr())
    plist = [os.path.join(outputdir, xfile)]
    outfile = os.path.join(outputdir, "results-j5evenodd-pareto-%s.oa" % adata0.idstr())
    if oapackage.checkFiles(outfile):
        continue
    parseParetoList(
        outfile, plist, verbose=0, afmode=oapackage.ATEXT, nrows=N, ncols=kk, paretomethod=paretomethod, cache=cache
    )
    npar = oapackage.nArrays(outfile)
    print("   Pareto arrays for column %d: %d arrays" % (kk, npar))


gjobs, jobs = jobStatus(alljobs)

# %%

joblistg1 = gatherResults([], outputdir, splitdata, adata=adata, dozip=False, verbose=2)
# [j.runjob() for j in joblistg1]

if 1:
    # j=gatherResults([], outputdir, splitdata, adata=adata, dozip=False, verbose=2)[0]
    j = joblistg1[0]
    j.ncores = 8
    print("-- top level gather command:")
    print(j)
    j.canrun(verbose=1)
    outfile, substr = createJobScript(j, queue="q72h", scriptdir=scriptdir)


# joblistg1[0].execute=True
# joblistg1[0].runjob()

alljobs += joblistg1

gjobs, jobs = jobStatus(alljobs)
print("----")
for j in gjobs:
    print("job: %s" % j)
print("----")
# STOP

# %%
# configfile='oaconfig%d.txt' % klevel

if args.statistics:
    (ntotal, nfiles) = evenoddAnalysePartialRun(outputdir, adata, splitdata)

print("Case: %s: script complete" % adata)


# %%
# TODO:
# 1. make sure j45 can be used everywhere (since j45 ordering can deviate from original ordering??!?!?!?)
# 4. incorporate pareto_develop2 (for web page generation)

if args.statistics:
    evenoddAnalyseRun(outputdir, adata, splitdata)


if 0:
    # compress files
    ll = findfiles(os.path.join(outputdir, edir), ".*oa$")
    if len(ll) > 0:
        cmd = "cd %s; gzip -7 -q *.oa" % os.path.join(outputdir, edir)
        os.system(cmd)

    if oahelper.checkFiles(os.path.join(outputdir, edir, mdfile), cache=cache) and 1:
        if verbose >= 2:
            print("md5 check already there")
    else:
        if verbose:
            print("creating md5 sums: %s" % os.path.join(outputdir, edir))
        researchOA.arrayfilesMD5(os.path.join(outputdir, edir), os.path.join(outputdir, edir, mdfile), verbose=1)


# %%

print("manual: oaextendsingle -l 1 -f B --maxk 8 -m 5")

# runcomplete, cmd = generateLevelTwo(ii, jj, splitdata, allcmdlogfile=allcmdlogfile)


# %%
if 0:

    lvls = [0]
    checkLevel(lvls, splitdata, adata, outputdir, verbose=1)


# %% Load packages


# gjobs, jobs = jobStatus(alljobs)


if args.queue is not None:
    queue = args.queue
else:
    queue = "q1h"
    queue = "q24h"
    queue = "q72h"
    queue = None
    # queue='q7d'


makeJobList(scriptdir, jobs)


runable_jobs = [j for j in alljobs if not j.complete() and j.canrun()]

if args.execute:
    print(f"## executing {len(runable_jobs)} jobs")
    job.execute = True
    [j.runjob() for j in runable_jobs]
    job_results = [j.runjob() for j in jobs]

if args.X == 2:
    if len(runable_jobs) > 0:
        # jx=alljobs[-1]
        print("last job: ")
        jx = runable_jobs[-1]
        print(jx)
        print(jx.checkfiles)

if args.X == 1:
    runable_jobs = [j for j in alljobs if not j.complete()]
    print("---- all jobs:")
    for jx in runable_jobs:
        jx.analyse(verbose=1)
        break

# %%

if 0:
    infile = join("/home/eendebakpt/misc/oa/oacode/build", "split-split-0.oa")
    outfile = join("/home/eendebakpt/misc/oa/oacode/build", "dummy.oa")
    # researchOA.parseParetoFile(join('/home/eendebakpt/misc/oa/oacode/build', 'split-split-0.oa'), outfile, oapackage.ABINARY)
    oapackage.calculateParetoEvenOdd([infile], outfile, 1, oapackage.ABINARY)
