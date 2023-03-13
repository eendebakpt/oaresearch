"""
Created on Fri Sep 21 16:14:21 2012

@author: eendebakpt
"""


# %% Load necessary packages
import copy
import itertools
import os
import pickle
import platform
import re
import shutil
import sys
import time
from os.path import basename, join
from typing import List, Union

import numpy as np
from oapackage.oahelper import safemax, series2htmlstr

from oaresearch.job_utils import job
from oaresearch.misc_utils import gzipOA

try:
    import matplotlib.pyplot as plt
except:
    # print('%s: matplotlib not available, not all functions will work...' %
    # __name__) # sys.modules[__name__]
    pass

import oapackage
import oapackage as oalib
import oapackage.markup as markup
import oapackage.oahelper as oahelper
import oapackage.scanf as scanf
from oapackage.markup import oneliner as e

from oaresearch.pythondevelop import ABhelper
from oaresearch.pythondevelop.ABhelper import *
from oaresearch.pythondevelop.ABhelper import checkFiles, checkOAfile

# %% Make nice plots
# http://blog.olgabotvinnik.com/prettyplotlib/


def niceplot(ax, fig=None, despine=True, verbose=0, legend=None, almost_black="#222222"):
    """Create a good looking plot

    The code:
        - removes spines
        - makes legend and spines lighter
        - makes legend lighter

    """

    # Remove top and right axes lines ("spines")
    if verbose:
        print("niceplot: remove spines")

    spines_to_keep = ["bottom", "left"]
    if despine:
        spines_to_remove = ["top", "right"]
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
    else:
        spines_to_keep += ["top", "right"]
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    if verbose:
        print("niceplot: reduce spine intensity")

    for spine in spines_to_keep:
        ax.spines[spine].set_linewidth(0.5)
        ax.spines[spine].set_color(almost_black)

    ax.tick_params(axis="both", direction="out")

    if legend is not None:
        if verbose:
            print("niceplot: adjust legend")

        # Remove the line around the legend box, and instead fill it with a light grey
        # Also only use one point for the scatterplot legend because the user will
        # get the idea after just one, they don't need three.
        light_grey = np.array([float(248) / float(255)] * 3)
        rect = legend.get_frame()
        rect.set_facecolor(light_grey)
        rect.set_linewidth(0.0)

        # Change the legend label colors to almost black, too
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)

    if fig is not None:
        fig.tight_layout(pad=1.0)

    plt.draw()
    plt.show()


# %% Document markup


def indexstart(ii=None):
    """Return starting index in papers (e.g. 0 or 1-based)"""
    if ii is None:
        return 1
    else:
        return ii + 1


def indexend(k):
    """Return end index in papers (e.g. 0 or 1-based)"""
    return k


# %% Math


def binom(n, k):
    """Return binomial coefficient"""
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


# %% File utilities


def findfiles(p, patt=None):
    """Get a list of files"""
    lst = os.listdir(p)
    if patt is not None:
        rr = re.compile(patt)
        lst = [l for l in lst if re.match(rr, l)]
    return lst


def finddirectories(p, patt=None):
    """Get a list of directories"""
    lst = os.listdir(p)
    if patt is not None:
        rr = re.compile(patt)
        lst = [l for l in lst if re.match(rr, l)]
    lst = [l for l in lst if os.path.isdir(os.path.join(p, l))]
    return lst


def splitBase(idx):
    """Return a named used when a big case is split into several pieces"""
    name = ""
    if len(idx) == 0:
        return "base"
    for ii, n in enumerate(idx):
        name += "sp%d-split-%d-" % (ii, n)
    name = name[0:-1]
    return name


def splitname(idx, oa=".oa"):
    """Return a named used when a big case is split into several pieces"""
    name = ""
    for ii, n in enumerate(idx):
        name += "sp%d-split-%d-" % (ii, n)
    name = name[0:-1]
    if oa is not None:
        name += oa
    return name


def splitfileTag(lvls, outputdir, reruntag):
    """Return name of zip file for splitted results"""
    edir0 = splitdir(lvls)
    splitfiletag0 = f"splitted-{splitTag(lvls)}-{reruntag}.txt"
    splitfiletag = os.path.join(outputdir, edir0, splitfiletag0)
    return splitfiletag0, splitfiletag


def splitTag(idx):
    """Return tag for given indices"""
    name = ""
    if len(idx) == 0:
        return "base"
    for ii, n in enumerate(idx):
        name += "%d." % n
    name = name[0:-1]
    return name


def splitdir(idx):
    """Return a subdirectory used when a big case is split into several pieces"""
    name = ""
    if len(idx) == 0:
        return ""
    for ii, n in enumerate(idx[0:]):
        name += "sp%d-split-%d/" % (ii, n)
    name = name[0:-1]
    return name


def mkdirc(d):
    """Make a directory if it does not already exist"""
    if not os.path.exists(d):
        try:
            os.mkdir(d)
        except:
            # catch for race conditions
            print("mkdirc: error with dir %s" % d)
    return d


def get_folder_size(path):
    """Return total size of direcotry in bytes"""
    TotalSize = 0
    for item in os.walk(path):
        for filex in item[2]:
            try:
                TotalSize = TotalSize + os.path.getsize(os.path.join(item[0], filex))
            except:
                print("error with file:  " + os.path.join(item[0], filex))
    return TotalSize


def analyseAB(afile, outputdir=None, Afinal=0, cache=1, verbose=0):
    """Calculate D-values and B-values for all arrays in a file"""
    adatan = oalib.arraydata_t(adata, kfinal)
    outfile = "finalselection-" + adatan.idstr() + "-A%.2f" % Afinal + ".oa"
    adatan = oalib.arraydata_t(adata, kfinal)

    anafiles = analyseFile(afile, method="full", verbose=1, cache=cache)
    anafile = anafiles[1]
    data = ABhelper.loadAnalysisFile(anafile, ncols=3)
    A = data[:, 1]
    B = data[:, 2]

    b = B
    b[b == 0] = np.inf

    # Final selection ##

    gidxmask = A >= Afinal
    bidxmask = A < Afinal
    gidx = gidxmask.nonzero()[0]
    bidxmask.nonzero()[0]

    safemax(A)
    safemin(b, np.inf)

    if outputdir is not None:
        join(outputdir, "tmp.oa")
        selectArrays(afile, outfile, gidx, verbose=1, cache=cache)


def splitFile(afile, ostr="split", nn=None, maxn=None, nwritemax=1000, method="a", verbose=1, cache=1, cmdlog=None):
    """Split an array file into several subfiles"""
    if nn is None:
        if maxn is None:
            nn = 10
        else:
            na = oahelper.nArrayFile(afile)
            nn = np.ceil(float(na) / maxn)
    cmd = "oasplit -i %s --nwritemax %d -n %d -o %s" % (afile, nwritemax, nn, ostr)
    if verbose >= 2:
        print("splitFile: cmd %s" % cmd)
    cfile = ostr + "-split-%d" % (0) + ".oa"
    if cmdlog is not None:
        cmdlog.write("# Split files\n" + cmd + "\n\n")
        cmdlog.flush()
    sys.stdout.flush()  # needed for tee command
    if not oahelper.checkFiles(cfile, cache):
        if verbose:
            print("splitting file %s in %d parts" % (basename(afile), nn))
        os.system(cmd)
        return ostr, nn
    else:
        if verbose >= 2:
            print("splitting file %s in %d parts (already done)" % (basename(afile), nn))
        return ostr, nn


def loadPickle(pfile):
    """Python 2/3 helper function"""
    try:
        rx = pickle.load(open(pfile, "rb"))
    except:
        rx = pickle.load(
            open(pfile, "rb"),
            encoding="latin",
        )

    return rx


def storeResults(rr, outfile, verbose=0):
    """Store results in Python pickle file"""
    rrs = copy.copy(rr)
    adata = rrs["adata"]
    if "fcopy" in rrs.keys():
        # print('storeResults')
        if "rr" in rrs["fcopy"]:
            rrs["fcopy"]["rr"].update({"adata": adata.idstrseriesfull()})
    if (type(adata)) == type(oalib.arraydata_t(1, 12, 1, 2)):
        if verbose:
            print("storeResults: replace adata structure")
        rrs["adata"] = adata.N, adata.ncols, adata.strength, adata.factor_levels()
    pickle.dump(rrs, open(outfile, "wb"))


def loadResults(infile):
    """Load results in Python pickle file"""
    rrs = loadPickle(infile)
    adv = rrs["adata"]
    ss = oalib.intVector(adv[3])
    rrs["adata"] = oalib.arraydata_t(ss, adv[0], adv[2], adv[1])
    return rrs


def loadbinfile(anafile, ncols=-1, verbose=0):
    data = np.fromfile(anafile, dtype=np.float64)
    nr = int(data.size / ncols)
    if data.size > 3:
        # we have a binary header
        # print('binary header')
        if data[0] == 30397995 and data[1] == 12224883:
            if verbose:
                print("loadbinfile: magic header")
            nr = data[2]
            if ncols != data[3] and ncols >= 0 and data[3] > 0:
                # unequal specification!
                raise NameError("loadbinfile %s: error nrows %d ncols %d data[3] %d" % (anafile, nr, ncols, data[3]))
            if ncols == -1:
                # take from data file
                ncols = data[3]

            data = data[4:]
    if ncols == -1:
        # empty array
        ncols = 0
    if verbose:
        print("data %s: %s: %d %s" % (anafile, data.shape, nr, ncols))
    data = data.reshape((int(nr), int(ncols)))
    return data


def arrayfilesMD5(xdir, outfile, verbose=1):
    """Create md5 sums for all array files in a directory"""
    ll = findfiles(xdir, ".*oa$")
    llg = findfiles(xdir, ".*oa.gz$")
    lst = ll + llg

    fid = open(outfile, "w")
    for name in lst:
        f = os.path.join(xdir, name)
        m = oalib.md5(f)
        fid.write(f"{m}  {name}\n")
    fid.close()


def arrayfilecompression(xdir, verbose=1):
    """Calculate compression factor"""
    ll = findfiles(xdir, ".*oa$")
    ll = ll + findfiles(xdir, ".*oa.gz$")
    tsize = 0
    ntotal = 0
    for name in ll:
        fname = os.path.join(xdir, name)
        sz = os.path.getsize(fname)
        nn = oalib.nArrays(fname)
        if verbose >= 2:
            print("file %s: size %d, n %d" % (name, sz, nn))
        tsize += sz
        ntotal += nn
    if ntotal == 0:
        c = np.NaN
    else:
        c = float(tsize) / ntotal
    if verbose:
        print("%d arrays, total size %.1f [GB]: %.1f bytes/array" % (ntotal, tsize / 1024.0**3, c))
    return c


def loadAnalysisFile(anafile, ncolshint=-1, verbose=1):
    """Load data from analysis file"""
    data = loadbinfile(anafile, ncols=ncolshint, verbose=verbose >= 2)
    data.shape[1]

    if verbose:
        print("loadAnalysisFile: loaded data %s" % str(data.shape))
    return data


def analyseFileCached(afile, method="full", verbose=1, cache=1):
    anafile = analyseFile(afile, method=method, verbose=1, cache=cache)
    if method == "full":
        data = loadAnalysisFile(anafile[1])
        Deff = data[:, 1]
        b = data[:, 2]
        rnk = data[:, 0]
        return data, Deff, b, rnk
    else:
        if verbose >= 2:
            print("analyseFileCached: file %s" % anafile[0])
        data = loadAnalysisFile(anafile[0])
        return data


def scriptCheck(cmd, lst):
    """Check whether all files in a list exist"""
    if not isinstance(lst, str):
        print("not implemented")
        return cmd
    cfile = lst

    scmd = "if [ ! -f %s ]; then\n" % cfile
    scmd += "  %s\n" % cmd
    scmd += "fi\n"
    return scmd


def dynamicExtendFile(afile, nextfile, kfinal, Afinal=0, verbose=1, dryrun=0, cache=1, logfile=None, cmdlog=None):
    """Extend arrayfile with dynamic filtering"""
    if kfinal is None:
        kfinal = -1  # fix this
    if logfile is None:
        cmd = "oaextendA -f D -d -v 2 -A %e -k %d %s -o %s" % (Afinal, kfinal, afile, nextfile)
    else:
        cmd = "oaextendA -f D -d -v 2 -A %e -k %d %s -o %s | tee %s" % (Afinal, kfinal, afile, nextfile, logfile)

    if verbose >= 2:
        print(cmd)
    if cmdlog is not None:
        cmdlog.write("# Extend dynamic\n" + scriptCheck(cmd, nextfile) + "\n\n")
        cmdlog.flush()
    sys.stdout.flush()  # needed for tee command
    if dryrun:
        return cmd

    # pdb.set_trace()#
    if not checkArrayFile(nextfile, cache=cache):
        stat = os.system(cmd)
        if stat:
            print('dynamicExtendFile: command "%s" failed!' % cmd)
            raise
    else:
        if verbose >= 2:
            icmd = "oainfo %s" % nextfile
            os.system(icmd)
            time.sleep(0.01)

    return cmd


def estimateNarrays(datadir, adata, maxn=None, cache=1, verbose=1):
    kfinal = adata.ncols
    k = kfinal
    os.chdir(datadir)
    adata.writeConfigFile("oaconfig.txt")
    t = adata.strength
    kstart = t + 1
    for ii in range(t, k + 1):
        adatax = oalib.arraydata_t(adata, ii)
        adatax.writeConfigFile("oaconfig%d.txt" % ii)

    if maxn is None:
        maxn = [1000] * (adata.ncols + 1)

    cmdlog = os.path.join(datadir, "commandlogn.txt")
    cmdlogfid = open(cmdlog, "w")

    # First analysis
    extendInitial(datadir, adata, kstart, verbose=verbose >= 2, cache=1, cmdlog=cmdlogfid)

    # Estimate total number of arrays

    sys.stdout.flush()
    if verbose >= 2:
        print("\n########## Estimate total number of arrays ##############")

    if not os.path.exists(os.path.join(datadir, "part")):
        os.mkdir(os.path.join(datadir, "part"))
    k = kstart
    adata0 = oalib.arraydata_t(adata, k)
    outfile0 = os.path.join(datadir, ("result-") + adata0.idstr() + ".oa")
    outfile = os.path.join(datadir, "part", "xx%d-" % k + adata0.idstr() + ".oa")
    shutil.copyfile(outfile0, outfile)
    shutil.copyfile(os.path.join(datadir, "oaconfig.txt"), os.path.join(datadir, "part", "oaconfig.txt"))

    rsplit = [1] * (kfinal + 1)
    rsplit[-1] = 0
    naa = [-1] * (kfinal + 1)
    nexpected = [-1] * (kfinal + 1)
    dt = [-1] * (kfinal + 1)
    for kk in range(kstart, kfinal):
        # Split
        # print('kk %d ' % kk)
        k = kk
        adata0 = oalib.arraydata_t(adata, k)
        adatan = oalib.arraydata_t(adata, k + 1)
        afile = os.path.join(datadir, "part", "xx%d-" % k + adata0.idstr() + ".oa")
        os.chdir(os.path.join(datadir, "part"))

        ostr = "xx%d" % (k)
        splitfile = ostr + "-split-0" + ".oa"
        na = ABhelper.nArrayFile(afile)
        naa[kk] = na
        nn = int(np.ceil(float(na) / maxn[k]))
        nexpected[kk] = int(np.ceil(float(na) / nn))

        ww = ABhelper.nArrayFile(splitfile)
        if not ww == nexpected[kk]:
            cache = 0
            if verbose:
                print(
                    "turning off cache! split factor %d splitfile %d nexpected %d na %d maxn[k] %d"
                    % (nn, ww, nexpected[kk], na, maxn[k])
                )
        if verbose >= 2:
            print(
                "  caching: split factor %d splitfile %d nexpected %d na %d maxn[k] %d"
                % (nn, ww, nexpected[kk], na, maxn[k])
            )

        (ostr, nn) = splitFile(
            afile, ostr=ostr, maxn=maxn[k], nwritemax=2, method="a", verbose=1, cache=cache, cmdlog=cmdlogfid
        )
        rsplit[kk] = nn
        ostrn = "xx%d" % (k + 1)
        nextfile = ostrn + "-" + adatan.idstr() + ".oa"

        # Extend
        # FIXME: what if A=0?

        logfile = nextfile.replace(".oa", "-log.txt")
        dynamicExtendFile(
            splitfile, nextfile, kfinal=kfinal, Afinal=0, verbose=1, cache=cache, logfile=logfile, cmdlog=cmdlogfid
        )
        dt[kk] = oahelper.parseProcessingTime(logfile)

        if verbose >= 2:
            print(
                "kk %d: %d arrays, split factor %d, %d arrays extended to %d arrays"
                % (kk, na, nn, nexpected[kk], nArrayFile(nextfile))
            )
        nlast = nArrayFile(nextfile)
        # print('afile %s, nextfile %s' % (afile, nextfile))

    naa[kfinal] = nArrayFile(nextfile)
    sfac = np.array(rsplit).prod()
    if verbose >= 2:
        print("estimated number of arrays %d = %.3e (split %d)" % (sfac * nlast, sfac * nlast, sfac))

    rr = dict(
        {
            "datadir": datadir,
            "rsplit": rsplit,
            "nexpected": nexpected,
            "narrays": naa,
            "adata": adata,
            "dt": dt,
            "kstart": kstart,
        }
    )
    return rr


def showNestimate(rr, verbose=1):
    adata = rr["adata"]
    rsplit = rr["rsplit"]
    kstart = rr["kstart"]
    print(adata.showstr())

    for ii in range(kstart, adata.ncols + 1):
        sfac = np.array(rsplit[kstart:ii]).prod()
        na = rr["narrays"][ii]
        dt = rr["dt"][ii]
        nstart = rr["nexpected"][ii]
        ne = sfac * na
        print(
            "%d columns: %d arrays, split factor %d (step %d), extending %d arrays, estimated number %.3e = %d"
            % (ii, na, sfac, rsplit[ii], nstart, ne, ne)
        )
        if verbose >= 2:
            print("   %d columns: processing in this step %.1f [s] = %.1f [m]" % (ii, dt, dt / 60.0))


def case2htmlstr2(ad, html=1):
    print("replace by series2htmlstr")
    s = list(ad.factor_levels())
    p = -1
    n = 0
    aa = []
    bb = []
    while n < len(s):
        if s[n] != p:
            p = s[n]
            aa.append(p)
            bb.append(1)
        else:
            bb[-1] += 1
        n = n + 1
    hstr = "OA(%d; %d; " % (ad.N, ad.strength)
    for ii in range(0, len(aa)):
        hstr += "%d<sup>%s</sup>" % (aa[ii], str(bb[ii]))
    hstr += ")"
    return hstr


def parseParetoFile(infile, outfile, verbose=1, afmode=oalib.ABINARY, nrows=None, ncols=None):
    """Calculate the Pareto optimal arrays from an array file

    Pareto optimality is calculated according to (rank; A3,A4; F4)
    """

    if nrows is None:
        nrows = 0
    if ncols is None:
        ncols = 0
    oapackage.calculateParetoEvenOdd([infile], outfile, verbose, afmode, nrows, ncols)

    return

    raise "old implementation"
    al = oalib.readarrayfile(infile)
    pset = oalib.parsePareto(al, 1)
    if verbose:
        pset.show(2)

    idx = pset.allindices()
    idx2 = oalib.intVector(idx)
    pal = oalib.arraylist_t()
    oalib.selectArrays(al, idx2, pal)
    if nrows is None or ncols is None:
        oalib.writearrayfile(outfile, pal, afmode)
    else:
        oalib.writearrayfile(outfile, pal, afmode, nrows, ncols)


def parseParetoList(
    outfile, lst, verbose=1, afmode=oalib.ATEXT, nrows=None, ncols=None, paretomethod=0, needall=True, cache=0
):
    """Calculate the Pareto optimal arrays from a list of array files"""

    if oahelper.checkFiles(outfile, cache=cache):
        if verbose:
            print("parseParetoList: file %s already exists" % outfile)
        return

    alist = oalib.arraylist_t()
    for ii, name in enumerate(lst):
        pfile = name
        if verbose >= 2:
            print("reading %s" % pfile)
        if oalib.nArrays(pfile) == -1:
            # problem with file, abort
            if needall:
                if verbose:
                    print("parseParetoList: missing file %s, not generating results" % pfile)
                return
            else:
                if verbose:
                    print("parseParetoList: missing file %s..." % pfile)
                continue
        oalib.readarrayfile(pfile, alist)

    if verbose >= 1:
        print("parseParetoList: calculating optimal designs from %d arrays" % alist.size())
    pset = oapackage.parsePareto(alist, verbose=1, paretomethod=paretomethod)
    pset.show(2)

    idx = pset.allindices()
    pal = oalib.arraylist_t()
    idx2 = oalib.intVector(idx)
    oalib.selectArrays(alist, idx2, pal)
    if nrows is None:
        oalib.writearrayfile(outfile, pal, afmode)
    else:
        oalib.writearrayfile(outfile, pal, afmode, nrows, ncols)


def loadnumbersfile(cfile, num=2, verbose=1):
    """Load text file with numbers
    Format: each line in the text file has the format:
        k %d: [%d]

    """
    fid = open(cfile)
    aa = np.zeros((0, 3))
    fmt = "k %d:" + " %d" * num
    for ln in fid:
        if verbose >= 2:
            print(ln)
        r = scanf.sscanf(ln, fmt)
        aa = np.vstack((aa, np.array(r).reshape((1, 1 + num))))
    fid.close()
    return aa


def combineNumbers(nlist, outfile=None, num=2, verbose=1):
    ww = dict()
    for nm in nlist:
        w = loadnumbersfile(nm)
        for ii, ln in enumerate(w):
            # print(ii)
            k = int(ln[0])
            if k in list(ww.keys()):
                ww[k] = ww[k] + ln[1:]
            else:
                ww[k] = ln[1:]
    if outfile is not None:
        fid = open(outfile, "w")
        for k in np.sort(list(ww.keys())):
            # ls='%d %d' %  (ww[k][0], ww[k][1])
            ls = " ".join(["%d" % x for x in ww[k]])
            fid.write("k %d: %s\n" % (k, ls))
        fid.close()
    return ww


def writeNumbers(nfile, data, verbose=1):
    """Write a numbers file to disk"""
    with open(nfile, "w") as fid:
        for k in np.sort(list(data.keys())):
            ls = " ".join(["%d" % x for x in data[k]])
            fid.write("k %d: %s\n" % (k, ls))


def citation(paper, style="brief"):
    """Return citation in html format
    Args:
        paper (str): paper to be cited
        style (str): brief or full
    """
    if paper == "complete":
        return markup.oneliner.a(
            "Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays",
            href="http://dx.doi.org/10.1002/jcd.20236",
        )
    elif paper == "conference" or paper == "cisomorphism":
        return "<em>A Classification Criterion for Definitive Screening Designs</em> (in preparation)"
        return markup.oneliner.a("A Classification Criterion for Definitive Screening Designs", href="...")
    elif paper == "conference enumeration" or paper == "cenumeration":
        return "<em>Enumeration and Classification of Definitive Screening Designs</em> (in preparation)"
        return markup.oneliner.a("Enumeration and Classification of Definitive Screening Designs", href="...")
    else:
        raise Exception("paper not known")


def parseParetoClass(adata, pfile, dstr=None, targetdir="", subfilebase="test", prevnextstr="", verbose=1):
    """Generate html page with data about arrays"""
    if verbose:
        print("parsing: %s" % pfile)

    strength = adata.strength
    xstrplain = oapackage.oahelper.series2htmlstr(adata, html=0, case=1)
    xstr = oapackage.oahelper.series2htmlstr(adata, html=1, case=1)
    al = oalib.readarrayfile(pfile)

    pset = oalib.parsePareto(al, 1)
    pset.show(2)

    """ Generate page """

    page = markup.page()

    page.init(
        title="Class %s" % xstrplain,
        css=("./oastyle.css"),
        lang="en",
        htmlattrs=dict({"xmlns": "http://www.w3.org/1999/xhtml", "xml:lang": "en"}),
        header="<!-- Start of page -->",
        bodyattrs=dict({"style": "padding-left: 3px;"}),
        # doctype=markup.doctype.strict,
        doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
        metainfo=(
            {
                "text/html": "charset=utf-8",
                "keywords": "orthogonal arrays designs",
                "robots": "index, follow",
                "description": "Even-Odd arrays",
            }
        ),
        footer="<!-- End of page -->",
    )

    if len(al) == 0:
        Fval = ""
        Fstr = ""
    else:
        A = al[0]
        js = oalib.jstruct_t(A, 4)
        Fval = js.Fval(strength)
        Fstr = "F(" + ",".join(["%d" % x for x in Fval]) + ")"

    page.h1("Class %s " % xstr)
    oap = markup.oneliner.a("Orthogonal Array", href="../software.html")
    pstr = "This page contains information about even-odd and foldover arrays. "
    pstr += "The results have been generated with the %s package." % oap
    pstr += (
        " If you use these data, please cite the paper "
        + markup.oneliner.a(
            "Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays",
            href="http://dx.doi.org/10.1002/jcd.20236",
        )
        + "."
    )
    page.p(pstr)

    m = 1 + adata.ncols + adata.ncols * (adata.ncols - 1) / 2

    if dstr is not None:
        page.p(dstr)
    page.p("For each array we have 3 different optimization criteria. These are: ")
    page.ol()
    tmp = "The rank of the second order interaction matrix.\n"
    tmp += "The maximum rank is %d." % m
    page.li(tmp)
    page.li("The values of A3 and A4 of the GWLP, followed by %s." % Fstr)
    page.li("The values %s," % Fstr)
    page.ol.close()
    tmp = "We have selected all Pareto optimal arrays according to these criteria.\n"
    if pset.numberindices() == 1:
        tmp += "There is <b>%d</b> Pareto optimal design." % (pset.numberindices())
    else:
        tmp += "There are <b>%d</b> Pareto optimal designs in <b>%d</b> different classes." % (
            pset.numberindices(),
            pset.number(),
        )
    page.p(tmp)

    page.table()
    page.tr(style="font-weight: bold; border-bottom: solid 1px black;")
    page.td(["Rank", "A3, A4; " + Fstr, Fstr, "Indices"], style="padding-right:30px;")
    page.tr.close()

    for ii, ee in enumerate(pset.elements):
        if verbose >= 2:
            print("ii %d" % ii)
        # print(ee.value)
        idx = ee.indices[0]
        A = al[idx]
        xf = oalib.array2xf(A)
        rnk = xf.rank()
        gwlp = A.GWLP()
        js = oalib.jstruct_t(A, 4)
        FF = js.calculateF(strength)
        Fval = js.Fval(strength)

        Fstr = ",".join(["%d" % x for x in FF])
        xfile = subfilebase + "-pareto-%d" % ii + ".oa"
        alsub = oalib.arraylist_t()
        oalib.selectArrays(al, ee.indices, alsub)
        oalib.writearrayfile(os.path.join(targetdir, xfile), alsub)

        if len(gwlp) < 5:
            gstr = "%.4f, -" % (gwlp[3]) + "; " + Fstr
        else:
            gstr = f"{gwlp[3]:.4f}, {gwlp[4]:.4f}" + "; " + Fstr
        page.tr(style="")
        page.td(["%d" % rnk, gstr, Fstr], style="padding-right:30px;")
        istr = ",".join(["%d" % x for x in ee.indices])
        istrlnk = markup.oneliner.a(istr, href=xfile)
        page.td(istrlnk, style="padding-right:3px;")
        page.tr.close()

    page.table.close()

    page.br(clear="both")
    page.p(prevnextstr)

    if 1:
        # statistis for all arrays
        page.br(clear="both")
        page.h2("Statistics for all Pareto arrays")

        xfile = subfilebase + "-pareto-all" + ".oa"
        lst = oalib.arraylist_t()
        for a in al:
            lst.push_back(a)

        oalib.writearrayfile(os.path.join(targetdir, xfile), lst)
        istrlnk = markup.oneliner.a("Pareto arrays", href=xfile)

        page.p("The file with all %s." % istrlnk)

        rrstats = generateStatistics(pfile, nbest=10000, verbose=1, doprojef=0)
        px = markup.page()
        px = statistics2htmltable(px, rrstats, verbose=1, titlestr="")
        page.p(str(px))

    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)
    page.br(clear="both")
    page.p("Page generated on %s." % dstr)

    if verbose >= 2:
        print(page)

    subfile = subfilebase + ".html"
    fid = open(os.path.join(targetdir, subfile), "w")
    fid.write(str(page))
    fid.close()

    return subfile


def oaseries2html(N, t, s=2, k="a"):
    """Return html string for design specification"""
    cstr = "OA(%d; %d; %d<sup>%s</sup>)" % (N, t, s, k)
    # cstr='OA(%d; %d<sup>%s</sup>; %d)' % (N, s, k, t)
    return cstr


def caseDefaultDir(ad):
    """Return case design string"""
    fullstr = "design-" + ad.idstr() + "-t%d" % ad.strength
    return fullstr


def hasResultsFiles(case, basedatadir, verbose=0):
    ad = case["case"]
    xdir = os.path.join(basedatadir, case["casedir"])
    for ii in range(ad.strength + 1, ad.ncols):
        adata0 = oalib.arraydata_t(ad, ii)
        afilebase = "result-%s.oa" % adata0.idstr()
        afile = os.path.join(xdir, afilebase)
        af = oalib.arrayfile_t(afile, 0)
        if not af.isopen():
            # print('problem with case %s' % case['casedir'])
            return False
    return True


def data2results(
    data, rr, nin=None, nout=None, ngoodrank=None, fullcalccurrent=None, fullcalc=None, selectionmaxn=None
):
    """Convert analysis data to results structure"""
    print("do not use this function!!!!")
    raise


def makeArraySelection(afile, gmaidx, outputfile, verbose=1, dogz=None, cache=1, hack=0):
    """Make a selection of arrays"""
    ngma = gmaidx.size
    gmaoutfile = outputfile
    if hack:
        print("makeArraySelection: debugging code")
        return (0, 0, 0, 0)
    if not (outputfile.endswith(".oa") or outputfile.endswith(".oa.gz")):
        raise NameError("makeArraySelection: filename %s" % outputfile)
    if verbose:
        print("  makeArraySelection: afile: %s, outputfile: %s, %d indices " % (afile, gmaoutfile, ngma))
    if dogz is None:
        dogz = gmaidx.size > 100
    if dogz:
        afmode = oalib.ABINARY
    else:
        afmode = oalib.ATEXT

    selectArrays(afile, gmaoutfile, gmaidx, afmode=afmode, verbose=1, cache=cache)
    if dogz:
        os.system("gzip -f %s" % gmaoutfile)
        outputfilefinal = outputfile + ".gz"
    else:
        outputfilefinal = outputfile
    if verbose >= 2:
        print("   makeArraySelection: done (dogz %d, afmode %d)" % (dogz, afmode))
    pxx = "arrays"
    if ngma == 1:
        pxx = "array"
    return (outputfilefinal, pxx, ngma, dogz)


def arraytxt(na):
    if na == 1:
        return "array"
    else:
        return "arrays"


def analyseFile(afile: str, method: str = "a", verbose=1, cache=1, cmdlog=None) -> Union[str, None, List]:
    """Analyse array file

    Args:
        afile: Array file to analyse
    Returns:
        Name of generated data file
    """
    nm = checkOAfile(afile)
    if nm is None:
        print("analyseFile: error: file %s does not exist" % afile)
        return None

    afile = nm
    if afile.endswith(".oa.gz"):
        repstr = ".oa.gz"
    else:
        repstr = ".oa"

    astr = afile.replace(repstr, "-ana")
    if verbose >= 2:
        print("analyseFile: method %s" % method)
    if method == "a":
        anafile: Union[str, List[str]] = afile.replace(repstr, "-ana-Dvalues.bin")
        cmd = f'echo "analyse file {afile}"; nice oaanalyse -j -1 -A -a {astr} {afile}'
    else:
        if verbose >= 2:
            print("analyseFile: method %s" % method)
        if method == "ab":
            anafile = [afile.replace(repstr, "-ana-Dvalues.bin"), afile.replace(repstr, "-ana-rankvalues.bin")]
            cmd = f'echo "analyse file {afile}"; nice oaanalyse -r -j -1 -A -a {astr} {afile}'
        elif method == "Deff":
            anafile = [afile.replace(repstr, "-ana-Dvalues.bin")]
            cmd = f'echo "analyse file {afile}"; nice oaanalyse -A -j -1 -a {astr} {afile}'
        elif method == "rankvalues":
            anafile = [afile.replace(repstr, "-ana-rankvalues.bin")]
            cmd = f'echo "analyse file {afile}"; nice oaanalyse -r -j -1 -a {astr} {afile}'
        elif method == "gma" or method == "gwlp":
            anafile = [afile.replace(repstr, "-ana-gwlpfull.bin")]
            cmd = f'echo "analyse file {afile}"; nice oaanalyse  -j -1 --gma -a "{astr}" "{afile}"'
        else:
            # full calculation
            anafile = [
                afile.replace(repstr, "-ana-Dvalues.bin"),
                afile.replace(repstr, "-ana-rankvalues.bin"),
                afile.replace(repstr, "-ana-gwlpfull.bin"),
            ]
            cmd = f'echo "analyse file {afile}"; nice oaanalyse --gma -r -j -1 -A -a {astr} {afile}'

    if verbose >= 2:
        print(cmd)
    # fid.write('# Analyse results\n' +  scriptCheck(cmd, anafile) + '\n\n')
    # fid.flush()
    if cmdlog:
        cmdlog.write("# Analyse results\n" + scriptCheck(cmd, anafile) + "\n\n")
        cmdlog.flush()
    sys.stdout.flush()  # needed for tee command
    if not checkFiles(anafile, cache=cache):
        if verbose >= 2:
            print(cmd)
        os.system(cmd)
    else:
        if verbose >= 2:
            print(f"file {anafile} already exists")
    return anafile


# %%


def array2html(
    X, header=1, tablestyle="border-collapse: collapse;", trclass="", tdstyle="", trstyle="", thstyle="", comment=None
):
    """Convert Numpy array to HTML table

    Arguments
    ---------
        X : numpy array
            array to be converted
        header : integer
            use header or not
    Returns
    -------
        page : markup html object
            generated table in HTML

    """
    page = markup.page()
    page.add("<!-- Created by array2html -->\n")
    if comment is not None:
        page.add("<!-- %s-->\n" % comment)

    page.table(style=tablestyle)
    offset = 0
    nc = X.shape[1]
    nr = X.shape[0]

    if isinstance(trstyle, str):
        trstyle = [trstyle] * nr
    if isinstance(trclass, str):
        trclass = [trclass] * nr

    ri = 0
    if header:
        page.tr(style="font-weight: bold; border-bottom: solid 1px black;" + trstyle[ri], class_=trclass[ri])
        ri = ri + 1
        for ii in range(nc):
            if isinstance(X[offset, ii], tuple):
                print("array2html: tuple instance")
                page.th(X[offset, ii][0], style=thstyle + X[offset, ii][1])
            else:
                page.th(X[offset, ii], style=thstyle)
        page.tr.close()
        offset = offset + 1

    nr = X.shape[0] - offset
    for r in range(nr):
        page.tr(style=trstyle[ri], _class=trclass[ri])
        for ii in range(nc):
            if isinstance(X[offset, ii], tuple):
                page.td(X[offset, ii][0], style=tdstyle + X[offset, ii][1])
            else:
                page.td(X[offset, ii], style=tdstyle)

        page.tr.close()
        offset = offset + 1
        ri = ri + 1
    page.table.close()
    return page


# %%


def loadGMAdata(afile, cache=1, verbose=0):
    af = oalib.arrayfile_t(afile)
    if af.isopen() == 0:
        del af
        return dict({"n": None, "indices": None, "gwp": None})
    afiles = analyseFile(afile, method="gma", verbose=verbose, cache=cache)
    # afiles[0]
    gmadata = loadAnalysisFile(afiles[0], ncolshint=1 + af.ncols, verbose=verbose)
    return gmadata


def analyseGMA(afile, cache=1, verbose=0):
    """Analyse an array file on GWLP values"""
    af = oalib.arrayfile_t(afile)
    if af.isopen() == 0:
        del af
        return dict({"n": None, "indices": None, "gwp": None})

    # FIXME: afile can be zipped, function cannot handle this...
    afiles = analyseFile(afile, method="gma", verbose=verbose, cache=cache)
    # afiles[0]
    gmadata = loadAnalysisFile(afiles[0], ncolshint=1 + af.ncols, verbose=verbose)
    if verbose >= 3:
        print(gmadata)

    # multiply by N and round to make numerically stable
    sind = sortrows(np.round(gmadata * af.nrows))
    if sind.size == 0:
        return dict({"n": 0, "indices": None, "gwp": None})

    vv = gmadata[sind, :]
    w = np.sum(np.abs(np.diff(np.round(vv * af.nrows), axis=0)), 1)
    w = np.hstack((w, [int(1)])).astype(int)
    nn = np.argmax(w != 0) + 1
    gmaidx = sind[0:nn]

    ii = sind[0]

    r = dict({"n": gmadata.shape[0], "indices": gmaidx, "gwp": gmadata[ii, :]})
    if gmadata.shape[0] < 100:
        # return all GWLP values
        r["fulldata"] = gmadata
    return r


def extendInitial(datadir, adata, kstart, verbose=1, cache=1, cmdlog=None):
    """Extend case to initial arrays"""

    # make sure root is present as well
    adata00 = oalib.arraydata_t(adata, adata.strength)
    outfile = os.path.join(datadir, "result-" + adata00.idstr() + ".oa")
    al = oalib.array_link(adata.N, adata.strength, -1)
    al.create_root(adata)
    oalib.writearrayfile(outfile, al)

    adata0 = oalib.arraydata_t(adata, kstart)
    cmd = "oaextendsingle -c oaconfig%d.txt -l 2 -f B -o result | tee log0.txt " % kstart
    if verbose:
        print("#extendInitial Initial calculation")
    if cmdlog is not None:
        cmdlog.write("# Initial calculation\n" + cmd + "\n\n")
        cmdlog.flush()
    cfile = "result-" + adata0.idstr() + ".oa"
    if verbose >= 2:
        print("cmd: %s" % cmd)
    sys.stdout.flush()  # needed for tee command
    if not checkFiles(cfile, cache):
        os.system(cmd)
    else:
        if verbose:
            print("#extendInitial: results already generated")

    adata0 = oalib.arraydata_t(adata, kstart)
    afile0 = os.path.join(datadir, "result-" + adata0.idstr() + ".oa")
    nsols = oahelper.nArrayFile(afile0)
    if verbose:
        print("  initial file: %d solutions" % nsols)
    return afile0, nsols


def projDeff(al, kp, verbose=1):
    """Calculate projection efficiencies of 2nd order models"""
    k = al.n_columns
    p = range(0, k)
    nc = oahelper.choose(k, kp)
    d = np.zeros(nc)
    pA = np.zeros(nc)
    for ii, c in enumerate(itertools.combinations(p, kp)):
        iv = oalib.intVector(c)
        als = al.selectColumns(iv)
        if verbose >= 2:
            print(c)
        d[ii] = als.Defficiency()
        pA[ii] = als.VIFefficiency()
    return (d, al.Defficiency(), pA)


# FIXME: make parsing structure for partial cases!


def generatePICtable(afile, kmin=3, kmax=5, verbose=1, aidx=None):
    alist = oalib.readarrayfile(os.path.join(afile))

    if aidx is None:
        aidx = range(0, len(alist))

    if verbose:
        print("generatePICtable: from array file %s" % afile)
    page = markup.page()
    page.add("\n")
    page.add("<!-- Created by generatePICtable -->\n")
    page.table(style=" border-collapse: collapse;")
    page.tr(style="font-weight: bold; border-bottom: solid 1px black;")
    page.th("Index", style="padding-right:30px; ")
    # page.th.close()
    page.th("D-efficiency", style="padding-right:10px; ")
    # page.th.close()
    for kp in range(kmin, kmax + 1):
        page.th("Mean %d-proj D-eff" % kp, style="padding-right:10px; ")
        page.th.close()
    page.tr.close()
    for ii, jj in enumerate(aidx):
        al = alist[jj]
        if verbose >= 2:
            print("generatePICtable: array %d->%d" % (ii, jj))
        page.add("   ")
        page.tr()
        page.td("%d" % jj)
        Deff = al.Defficiency()
        page.td("%.4f" % Deff)
        for kp in range(kmin, kmax + 1):
            (Dproj, Deff, Aprof) = projDeff(al, kp=kp)
            page.td("%.4f" % np.mean(Dproj))
        page.tr.close()
    page.table.close()
    page.add("\n")
    return page


def statistics2htmltable(page, rrstats, verbose=1, titlestr=None):
    if rrstats is not None:
        data = rrstats["data"]
        if "gwlpdata" in list(rrstats.keys()):
            gwlp = rrstats["gwlpdata"]
            # gwpl=rrstats['gma']['fulldata']
        else:
            gwlp = None
        na = data.shape[0]
        if verbose:
            print("statistics2htmltable:  generating full statistics table".format())
        if titlestr is None:
            if rrstats["mode"] == "all":
                page.h2("Table of arrays")
            else:
                page.h2("Table of arrays (selection)")
        else:
            if titlestr == "":
                if verbose:
                    print("statistics2htmltable: no title")
            else:
                page.h2(titlestr)

        # page.p()
        page.add("<!-- Created by statistics2htmltable -->\n")
        page.table(style=" border-collapse: collapse;")
        page.tr(style="font-weight: bold; border-bottom: solid 1px black;")
        page.th("Array", style="padding-right:30px; ")
        page.th(("Rank (2nd order)", "D-efficiency", "Average VIF"), style="padding-right:14px;")
        page.th(("E-efficiency"), style="padding-right:14px;")
        page.th(("GWLP"), style="padding-right:14px;")
        page.tr.close()
        iii = rrstats["aidx"]
        for ii in range(0, na):
            page.add("  ")
            page.tr(style="font-weight: normal;")
            page.td("%d" % iii[ii], style="padding-right:10px;")
            page.td("%d" % data[ii, 0], style="padding-right:1px;")
            xft = "%.5f"
            page.td(xft % data[ii, 1], style="padding-right:1px;")  # D-eff
            if data[ii, 2] == 0:
                page.td("Inf", style="padding-right:1px;")
            else:
                page.td(xft % data[ii, 2], style="padding-right:1px;")
            page.td(xft % data[ii, 3], style="padding-right:1px;")
            if gwlp is None:
                gstr = "-"
            else:
                gstr = oahelper.gwlp2str(gwlp[ii, :])
            page.td(e.small(gstr), style="padding-right:1px;")
            page.tr.close()
        page.table.close()
        # page.p.close()
        return page


def formatAtag(txt, lnk, afile):
    """Create HTML link to .oa file

    Args:
        txt (str)
        lnk (str)
        afile (str): .oa file on disk
    Returns:
        str: generated html

    """
    if oahelper.oaIsBinary(afile):
        ss = e.a(txt + e.small(" (binary)"), href=lnk, class_="binarylink")
    else:
        ss = e.a(txt, href=lnk, class_="normal")
    return ss


def abSubpage(rr, htmldir, verbose=1, makeheader=True, cache=1):
    """Generate html page for D-efficiency analysis"""
    ad = rr["adata"]
    idstr = "class-%s-t%d" % (ad.idstr(), ad.strength)
    if makeheader:
        pext = "html"
    else:
        pext = "html.template"
    if rr["mode"] == "special":
        subfile = "classdata-special-%s-t%d.%s" % (ad.idstr(), ad.strength, pext)
        htmlsubdir0 = "classdata-special-%s-t%d" % (ad.idstr(), ad.strength)
    else:
        subfile = "classdata-%s-t%d.%s" % (ad.idstr(), ad.strength, pext)
        htmlsubdir0 = "abdata-%s-t%d" % (ad.idstr(), ad.strength)
    rr["htmlsubdir0"] = htmlsubdir0
    xstr = series2htmlstr(ad, case=1)
    xstrplain = series2htmlstr(ad, html=0, case=1)

    if verbose:
        print("abSubpage: start of generation")

    htmlsubdir = os.path.join(htmldir, "tpages", htmlsubdir0)
    if not os.path.exists(htmlsubdir):
        os.mkdir(htmlsubdir)

    if not os.path.exists(htmlsubdir):
        os.mkdir(htmlsubdir)

    subfilef = os.path.join(htmldir, "tpages", htmlsubdir0, subfile)
    if checkFiles(subfilef, cache):
        if verbose:
            print("abSubpage: %s already exists" % subfile)
        return subfile
    else:
        if verbose:
            print("abSubpage: data: %s" % (xstr))
    datadir = rr["datadir"]

    page = markup.page()
    if makeheader:
        page.init(
            title="Class %s" % xstrplain,
            css=("../oastyle.css"),
            lang="en",
            htmlattrs=dict({"xmlns": "http://www.w3.org/1999/xhtml", "xml:lang": "en"}),
            header="<!-- Start of page -->",
            bodyattrs=dict({"style": "padding-left: 3px;"}),
            # doctype=markup.doctype.strict,
            doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
            metainfo=(
                {
                    "text/html": "charset=utf-8",
                    "keywords": "orthogonal arrays designs",
                    "robots": "index, follow",
                    "description": "Complete Enumeration of Orthogornal Arrays",
                }
            ),
            footer="<!-- End of page -->",
        )

    page.h1("Class %s " % xstr)
    oap = e.a("Orthogonal Array", href="../../software.html")
    pstr = "This page contains information about D-efficiency and GWLP values of the arrays. "
    pstr += "The results have been generated with the %s package." % oap
    pstr += (
        " If you use these data, please cite the paper "
        + e.a(
            "Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays",
            href="http://dx.doi.org/10.1002/jcd.20236",
        )
        + "."
    )
    page.p(pstr)
    # 1. calculate values

    if rr["mode"] == "full":
        page.h2("Results")
    else:
        if rr["mode"] == "special":
            page.h2("Results (special calculation)")
        else:
            page.h2("Results (partial calculation)")

    page.table()
    page.tr(style="font-weight: bold;")
    page.td("Statistic", style="padding-right:30px;")
    page.td(("Results"), style="padding-right:8px;")
    page.tr.close()

    page.tr(style="")
    page.td("Number of non-isomorphic arrays", style="padding-right:30px;")
    page.td(str(rr["narrays"]), style="padding-right:8px;")
    page.tr.close()

    gstr = rr["gmastr"]
    page.tr(style="")
    page.td("Best GWLP", style="padding-right:22px;")
    page.td(gstr, style="padding-right:8px;")
    page.tr.close()

    # rstr=('%d' % rr['narrays'], '%.4f' % rr['maxA'], '%.4f' % rr['minB'],
    # xxx, gstr)

    afstr = ""
    bfstr = ""
    if rr["ABselection"] is not None:
        if verbose >= 2:
            print("abSubpage: make AB selection")
        ww = rr["ABselection"]
        shutil.copyfile(os.path.join(datadir, ww["Afile"]), os.path.join(htmlsubdir, ww["Afile"]))
        shutil.copyfile(os.path.join(datadir, ww["Bfile"]), os.path.join(htmlsubdir, ww["Bfile"]))
        # print('abSubpage: copy: %s' % ww['Bfile'] )

        txt = "[%d %s with highest D-efficiency value]" % (ww["nA"], arraytxt(ww["nA"]))
        lnk1 = formatAtag(txt, ww["Afile"], os.path.join(datadir, ww["Afile"]))
        # pdb.set_trace()

        afstr = e.small("&nbsp;&nbsp;" + lnk1)

        txt = "[%d %s with best average VIF value]" % (ww["nB"], arraytxt(ww["nB"]))
        lnk1 = formatAtag(txt, ww["Bfile"], os.path.join(datadir, ww["Bfile"]))
        # pdb.set_trace()
        # print('grrr')
        bfstr = e.small("&nbsp;&nbsp;" + lnk1)
        # bfstr=e.a(' %d arryays' % ww['nB'], href=ww['Bfile'])
        # pdb.set_trace()

    if "gmafile" in list(rr.keys()):
        if verbose >= 2:
            print("abSubpage: make gmafile")

        if rr["gmafile"] is not None:
            shutil.copyfile(os.path.join(datadir, rr["gmafile"]), os.path.join(htmlsubdir, rr["gmafile"]))

    def simpleRow(a, b):
        page.tr(style="")
        page.td(a, style="padding-right:30px;")
        page.td(b, style="padding-right:8px;")
        page.tr.close()

    if rr["mode"] == "full":
        simpleRow("D-efficiency", ("Best value: %.4f " % rr["maxDeff"] + afstr))
        simpleRow("Average VIF", ("Best value %.4f " % rr["minB"] + bfstr))
        simpleRow("E-efficiency", ("Best value %.4f " % rr["maxE"]))
    else:
        simpleRow("D-efficiency", ("Value: %s " % str(rr["maxDeff"]) + afstr))

    havearrayfile = 0
    if "arrayfile" in list(rr.keys()):
        if rr["arrayfile"] is not None:
            havearrayfile = 1

    if havearrayfile:
        iafile = rr["arrayfile"]
        outfile0 = idstr + "selection" + ".oa"

        tmpfile = os.path.join(rr["datadir"], rr["arrayfile"])
        # print('hack tmpfile %s ' % tmpfile );  na=-1 # hack
        na = nArrayFile(tmpfile)
        if verbose >= 2:
            print("abSubpage: make arrayfile: na %d" % na)
        # print('hack %s ' % os.path.join( rr['datadir'], rr['arrayfile']) );  na=-1 # hack
        # pdb.set_trace()

        if na < 20000 and na >= 0:
            if iafile.endswith(".gz"):
                outfile0 += ".gz"
            outfile0final = copyOAfile(
                os.path.join(rr["datadir"], iafile), htmlsubdir, outfile0, convert=100, zipfile=None, verbose=0, cache=0
            )
            # shutil.copyfile(os.path.join( rr['datadir'], rr['arrayfile']),
            # os.path.join( htmlsubdir,outfile0) )
            if rr["fullcalc"]:
                # pdb.set_trace()
                xxx = formatAtag("all arrays", outfile0final, os.path.join(htmlsubdir, outfile0final))

                # xxx=e.a('all arrays', href=outfile0final)
                rr["datafilestr"] = "all arrays"
            else:
                na = nArrayFile(os.path.join(htmlsubdir, outfile0final))
                xxx = formatAtag("%d array(s)" % na, outfile0final, os.path.join(htmlsubdir, outfile0final))
                # xxx=e.a('%d array(s)' % na, href=outfile0)
                rr["datafilestr"] = "%d array(s)" % na

            page.tr(style="")
            page.td("Data", style="padding-right:30px;")
            page.td(xxx, style="padding-right:8px;")
            page.tr.close()

            rr["datafilestrfile"] = outfile0final
        elif na > 0:
            # pdb.set_trace()
            outputfile = os.path.join(htmlsubdir, outfile0)
            if verbose >= 2:
                print(
                    "abSubpage: make arrayfile: selectBestArrays (arrayfile %s, outputfile %s)"
                    % (rr["arrayfile"], outputfile)
                )
            sd = selectBestArrays(rr["datadir"], rr["arrayfile"], outputfile, nmax=2000, cache=1, verbose=0)
            #            sd=dict('')
            #            sd['filename']='-'
            #            sd['dstr']='-'; xxx=''; na=-1

            if verbose >= 2:
                print("abSubPage: after selectBestArrays")

            rr["datafilestrfile"] = basename(sd["filename"])
            rr["datafilestr"] = "best " + sd["dstr"]
            rr["datafilestr"] = "" + sd["dstr"]
            na = sd["aax"][2]
            xxx = formatAtag("%d array(s)" % na, rr["datafilestrfile"], outputfile)

            #            xxx=e.a('%d array(s)' % na, href=rr['datafilestrfile'])

            page.tr(style="")
            page.td("Data", style="padding-right:30px;")
            page.td(xxx, style="padding-right:8px;")
            page.tr.close()

        else:
            if verbose:
                print("abSubpage: no datafile (na %d)" % na)
            rr["datafilestrfile"] = None
            rr["datafilestr"] = "-"

    if "outfileVIF" in list(rr.keys()) and rr["outfileVIF"] is not None:
        print("generating page for VIF arrayfile")

        if verbose >= 3:
            print("abSubpage: make VIF arrayfile")

        iafile = rr["outfileVIF"]
        outfile0 = idstr + "selectionVIF" + ".oa"
        # print('-----------\n abSubpage: make VIF arrayfile: %s ' % iafile)

        tmpfile = os.path.join(rr["datadir"], rr["outfileVIF"])
        na = nArrayFile(tmpfile)

        outfile0final = copyOAfile(tmpfile, htmlsubdir, outfile0, convert=100, zipfile=100, verbose=1, cache=0)

        xxx = formatAtag("%d array(s)" % na, basename(outfile0final), os.path.join(htmlsubdir, outfile0final))
        #        xxx=e.a('%d array(s)' % na, href=basename(outfile0final) )

        page.tr(style="")
        page.td("Data (%d best VIF arrays)" % na, style="padding-right:30px;")
        page.td(xxx, style="padding-right:8px;")
        page.tr.close()
    else:
        if verbose >= 2:
            print("no VIF arrayfile ")

    if verbose >= 3:
        print("  abSubPage: at totaltime")

    if "totaltime" in list(rr.keys()):
        page.tr(style="")
        page.td("Processing time", style="padding-right:30px;")
        if isinstance(rr["totaltime"], bytes):
            page.td("%s" % rr["totaltime"], style="padding-right:8px;")
        else:
            if rr["totaltime"] < 60:
                page.td(" &lt; 1 minute ", style="padding-right:8px;")
            elif rr["totaltime"] < 3600:
                page.td("%.1f minutes (estimate)" % (float(rr["totaltime"]) / 60.0), style="padding-right:8px;")
            else:
                page.td("%.1f hours (estimate)" % (float(rr["totaltime"]) / 3600.0), style="padding-right:8px;")
        page.tr.close()

    page.table.close()

    if rr["abscatterplot"] is not None:
        subimage = rr["abscatterplot"]
        if verbose:
            print("  generating scatterplot image: %s" % subimage)
        figfile = rr["abscatterplot"]
        shutil.copy(os.path.join(datadir, figfile), os.path.join(htmlsubdir, figfile))
        if verbose >= 2:
            print(
                "  generating scatterplot image: copy %s to %s"
                % (os.path.join(datadir, figfile), os.path.join(htmlsubdir, figfile))
            )
        page.h2("Scatter plot")
        page.p()
        # page.img( src=htmlsubdir0 + '/' + subimage, style='max-width: 80%')
        page.img(src=subimage, style="max-width: 80%", alt="")
        page.p.close()

    statistics2htmltable(page, rr["fullstatistics"])

    if "projectiontable" in rr.keys():
        pp = rr["projectiontable"]
        page.add("<br/>\n")
        page.add(str(pp))

    # import datetime
    # today = datetime.date.today()
    # print today
    # print today.strftime('We are the %d, %h %Y')
    # today.strftime('%c')
    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)
    # print "Local current time :", localtime

    page.p("<br/>\n")  # clear='left')
    page.p("Page generated on %s." % dstr)
    # page.p.close()
    # subpage: raw data
    # page.pre(xx)

    pstr = str(page).replace(
        '<meta content="charset=utf-8" name="text/html" />',
        '<meta http-equiv="Content-Type" content="charset=utf-8" name="text/html" />',
    )

    print("writing to %s" % subfilef)
    fid = open(subfilef, "w")
    fid.write(pstr)
    fid.close()

    if verbose >= 3:
        print(f"abSubpage: {subfile} ({subfilef})")
    return subfile


def formatProccessingTime(ss, verbose: int = 1, estimate: bool = True, keep_seconds=False):
    """Format processing time to string

    Args:
        ss: Time in seconds or a string
    """
    if isinstance(ss, (str, bytes)):
        res = ss
    else:
        if ss < 0:
            res = "-1"
        elif ss < 60:
            if keep_seconds:
                res = "%.1f seconds" % (float(ss))
            else:
                res = "&lt; 1 minute"
        elif ss < 3600:
            res = "%.1f minutes" % (float(ss) / 60.0)
            if estimate:
                res += " (estimate)"
        else:
            res = "%.1f hours" % (float(ss) / 3600.0)
            if estimate:
                res += " (estimate)"
    return res


def N2kmax(N):
    """Calculate maximum k value for N"""
    a = 0.5
    b = 0.5
    c = 1 - N
    D = b**2 - 4 * a * c
    k = (-b + np.sqrt(D)) / (2 * a)
    k = int(np.floor(k))
    return k
    # 1+k+k(k-1)/2=N


def k2nmin(k, t=2, verbose=0):
    """Calculate minimum value of N for a given k and t such that a second order model can be estimated"""
    N = 1 + k + k * (k - 1) / 2
    Nmin = (2**t) * np.ceil(float(N) / 2**t)
    if verbose:
        print(N)
        print(Nmin)
    return int(Nmin)


def seriesSubpage(basedatadir, htmldir, case, dogma=0, verbose=1, cache=1, subtag="series", makepageheader=True):
    ad = case["case"]
    subfile = "series-%s-t%d.html" % (ad.idstr(), ad.strength)
    if checkFiles(os.path.join(htmldir, subtag, subfile), cache):
        if verbose:
            print("seriesSubpage: %s already exists" % subfile)
        return subfile, None

    numbersfile = "numbers-%s-t%d.txt" % (ad.idstr(), ad.strength)
    nums = caseReadNumbers(case, basedatadir)
    xstr = series2htmlstr(ad)
    xdir = os.path.join(basedatadir, case["casedir"])
    xstrplain = series2htmlstr(ad, html=0)

    page = markup.page()
    if makepageheader:
        page.init(
            title="Series %s" % xstrplain,
            css=("oastyle.css"),
            lang="en",
            header="",
            htmlattrs=dict({"xmlns": "http://www.w3.org/1999/xhtml"}),
            # doctype=markup.doctype.strict,
            doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"   "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">',
            metainfo=(
                {
                    "keywords": "orthogonal arrays designs",
                    "robots": "index, follow",
                    "description": "Complete Enumeration of Orthogornal Arrays",
                }
            ),
            footer="<!-- End of page -->",
        )

    page.h1("Series %s " % xstr)
    oap = e.a("Orthogonal Array", href="http://www.pietereendebak.nl/oapage/software.html")
    oap = e.a("Orthogonal Array", href="../software.html")
    if dogma:
        page.p(
            "This page contains information about the GWLP (generalized word length pattern) of the GMA design in a class. The results have been generated with the %s package."
            % oap
        )
    else:
        page.p(
            "This page contains information tne number of arrays in the series. The results have been generated with the %s package."
            % oap
        )

    # 1. calculate values

    print("dir: %s, dogma: %d" % (xdir, dogma))
    page.h2("Number of arrays")
    page.table()
    page.tr(style="font-weight: bold;")
    page.td("Case", style="padding-right:20x;")
    if dogma:
        page.td(("Number of arrays", "Best GWLP"), style="padding-right:8px;")
    else:
        page.td(("Number of arrays"), style="padding-right:8px;")
    page.tr.close()

    kmax = min(ad.ncols, len(nums) + ad.strength)
    for ii in range(ad.strength + 1, kmax + 1):
        if verbose:
            print("seriesSubPage: ii %d" % ii)
        adata0 = oalib.arraydata_t(ad, ii)
        xstr = series2htmlstr(adata0, case=1)
        afilebase = "result-%s.oa" % adata0.idstr()
        gmafilebase = "gmaarray-t%d-%s-gma.oa" % (adata0.strength, adata0.idstr())

        if verbose >= 2:
            print("file %s" % afilebase)
        afile = os.path.join(xdir, afilebase)
        af = oalib.arrayfile_t(afile, 0)
        vv = af.isopen()
        nv = af.narrays
        af.closefile()
        dorow = dogma and (vv) and (nv < 8000)
        if dorow:
            if ABhelper.checkFiles(afile, cache):
                if verbose >= 2:
                    print("file exists")
            else:
                if verbose >= 2:
                    print("  analyse GMA")
                # continue
            # FIXME: afile can be zipped
            rr = analyseGMA(afile, cache=cache)
            na = rr["n"]
            gmastr = ABhelper.gma2str(rr["gwp"], t=None)
            gmaidx = rr["indices"]
            if na > 0:
                if verbose >= 2:
                    print("  analyse GMA: generate file!")
                    if ii == 5:
                        pdb.set_trace()
                gmaoutfile = os.path.join(htmldir, subtag, gmafilebase)
                dogz = gmaidx.size > 300 or (gmaidx.size > 80 and ad.N > 91)
                gmaoutfilefinal, pxx, ngma, tmp = makeArraySelection(
                    afile, gmaidx, gmaoutfile, cache=cache, verbose=1, dogz=dogz
                )
                if dogz:
                    gmafilebasex = gmafilebase + ".gz"
                else:
                    gmafilebasex = gmafilebase
                gmastr += " " + e.a("(%d %s)" % (ngma, pxx), href="%s" % gmafilebasex)
            else:
                gmastr += " "

            print("updated gma string")
        else:
            na = nums[ii - ad.strength - 1]
            gmastr = "-"
        page.tr()
        page.td(xstr, style="padding-right:14px;")
        if dogma:
            page.td((str(na), gmastr))
        else:
            page.td(str(na))
        page.tr.close()
        # page.p( xstr + ': ' + 'gma: %s' % gmastr)

    page.table.close()

    # subpage: make image
    if len(nums) > 5:
        if verbose:
            print("generating image")
        page.h2("Image")
        xl = list(range(ad.strength + 1, len(nums) + ad.strength + 1))
        pnums = np.array(nums, dtype=np.float)
        pnums[pnums < 0] = np.nan
        plt.figure(1)
        plt.clf()
        plt.plot(xl, pnums, ".b", markersize=20)
        plt.title(r"$%s$" % ad.latexstr(), fontsize=18)
        plt.ylabel("Number of arrays", fontsize=15)
        plt.xlabel("Number of columns", fontsize=15)
        plt.xticks(xl, xl)
        plt.ylim(-1, max(nums) + 1)
        plt.xlim(min(xl) - 1, max(xl) + 1)

        subimage = "series-%s-t%d.png" % (ad.idstr(), ad.strength)
        plt.savefig(os.path.join(htmldir, "series", subimage))

        page.img(src="" + subimage, style="max-width: 70%")
    #                    impage.div(style='float: right; max-width: 50%; min-width: 300px;')
    #           impage.img(src=imbase, width='')

    # subpage: raw data
    page.h2("Raw processing file")
    # fid=open(os.path.join(htmldir, 'files', numbersfile) )
    fid = open(os.path.join(xdir, numbersfile))
    xx = fid.read()
    fid.close()
    page.pre(xx)

    return subfile, page


def selectionArraysScore(afile, outfile, vals, verbose=1, nmax=100, order="ascending", cache=1):
    """Select arrays from a file base on score"""
    vals = vals.flatten()
    idx = np.argsort(vals)

    if order == "descending":
        idx = idx[::-1]

    n = np.min((nmax, idx.size))
    gidx = idx[0:n]

    if n > 100:
        afmode = oalib.ABINARY
    else:
        afmode = oalib.ATEXT

    selectArrays(afile, outfile, gidx, verbose=verbose, cache=cache, afmode=afmode)

    return (afile, outfile, n)


def generateABcase(
    N,
    kfinal,
    kstart=3,
    kfull=3,
    Afinal=0,
    t=2,
    verbose=1,
    pfig=None,
    cache=1,
    selectfullrank=0,
    basedir="/home/eendebakpt/oatmp/final/",
    oadir="/home/eendebakpt/misc/oa/oacode/",
):
    """Helper function: extend a full series with D-efficiency selection"""

    datadir = basedir + "/oa%d-%d-t%d-small-fullrank%d" % (N, kfinal, t, selectfullrank)
    if not os.path.exists(datadir):
        os.mkdir(datadir)
    os.chdir(datadir)

    k = kfinal
    m = 1 + k + k * (k - 1) / 2
    s = oalib.intVector([2] * kfinal)
    adata = oalib.arraydata_t(s, N, t, k)
    if verbose:
        print("generateABcase: case: %s, k %d, m %d" % (str(adata), k, m))
    selectionmaxn = [1e9] * kfinal

    adata.writeConfigFile("oaconfig.txt")
    for ii in range(t, k + 1):
        adatax = oalib.arraydata_t(adata, ii)
        adatax.writeConfigFile("oaconfig%d.txt" % ii)

    cmdlog = os.path.join(datadir, "commandlog.txt")
    cmdlogfid = open(cmdlog, "w")

    # Create root array
    adata0 = oalib.arraydata_t(adata, t)
    al = oalib.array_link(N, t, -1)
    al.create_root(adata0)
    sols = oalib.arraylist_t()
    sols.append(al)
    outfile = ("result-") + adata0.idstr() + ".oa"
    oalib.writearrayfile(outfile, sols)
    cmdlogfid.write("# Create root array, write to %s" % outfile)

    # First analysis
    extendInitial(datadir, adata, kstart, verbose=1, cache=1, cmdlog=cmdlogfid)

    adata0 = oalib.arraydata_t(adata, kstart)
    afile0 = os.path.join(datadir, "result-" + adata0.idstr() + ".oa")
    if verbose:
        print("  initial file: %d solutions" % oalib.nArrays(afile0))

    outfile = ("result-") + adata0.idstr() + ".oa"
    anafile = analyseFile(outfile, method="full", verbose=1, cache=1)
    anafile = anafile[1]

    # init
    k = kstart
    adata0 = oalib.arraydata_t(adata, kstart)
    outfile0 = ("result-") + adata0.idstr() + ".oa"
    outfile = "resultdynamic-" + adata0.idstr() + ".oa"
    shutil.copyfile(outfile0, outfile)

    totaltime = 0

    for kk in range(kstart, kfull):
        # Initial input: show
        k = kk
        adata0 = oalib.arraydata_t(adata, k)
        afile = os.path.join(datadir, "resultdynamic-" + adata0.idstr() + ".oa")
        anafile = analyseFile(afile, verbose=1, cache=cache)
        # Calculate: k to k+1
        kn = k + 1
        adatan = oalib.arraydata_t(adata, kn)
        nextfile = os.path.join(datadir, "resultdynamic-" + adatan.idstr() + ".oa")
        logfile = nextfile.replace(".oa", "-log.txt")
        dynamicExtendFile(
            afile, nextfile, kfinal=kfinal, Afinal=Afinal, cmdlog=cmdlogfid, verbose=1, logfile=logfile, cache=cache
        )
        dt = parseProcessingTime(logfile)
        totaltime += dt
        print(f"processing: {dt:.0f} [s], {float(dt) / 3600:.1f} [h]")

    adata0 = oalib.arraydata_t(adata, kfull)
    infile = "resultdynamic-" + adata0.idstr() + ".oa"
    nextfile = "selecteddynamicresult-" + adata0.idstr() + ".oa"
    shutil.copyfile(infile, nextfile)

    # Prepare for selection rounds
    for kk in range(kfull, kfinal):
        if verbose:
            print("# selection round %d->%d" % (kk, kk + 1))

        k = kk
        adata0 = oalib.arraydata_t(adata, k)
        afile = os.path.join(datadir, "selecteddynamicresult-" + adata0.idstr() + ".oa")
        outfile = "selecteddynamic-" + adata0.idstr() + ".oa"
        m = int(1 + adata0.ncols + adata0.ncols * (adata0.ncols - 1) / 2)
        if checkFiles(outfile, cache=cache):
            if verbose >= 3:
                print("  skipping analysis")
            cmdlogfid.write("# Select best [x] arrays, write to file %s\n\n" % (outfile))
        else:

            if selectfullrank:
                anafiles = analyseFile(afile, method="full", verbose=1, cache=cache)
                data = loadAnalysisFile(anafiles[1])
                a = data[:, 1]
                r = data[:, 0]
                ngoodrank = (r == m).sum()
                ngoodrank2 = (a > 1e-10).sum()
                if ngoodrank != ngoodrank2:
                    print("warning: rank calculation has numerical stability issues?")
            else:
                anafiles = analyseFile(afile, verbose=1, method="full", cache=cache)
                data = loadAnalysisFile(anafiles[1])
                a = data[:, 1]
            drawAvalues(data, fig=100 + k)
            plotAthresholdsY(Afinal, kfinal, k)
            plt.title(r"Case $%s$: %d columns, selection" % (adata.latexstr(), k))

            if 1:
                Dvals = data[:, 1]
                Amax = oahelper.safemax(Dvals)
                Amin = oahelper.safemin(Dvals)
                #                rr[kk]['Amax']=Amax
                print("generateABcase: %d columns: Amax %.4f, Amin %4f", (kk, Amax, Amin))

            if data.size > 0:
                fraction = min(max(float(selectionmaxn[k]) / data.size, 0.0001), 1)
            else:
                fraction = 1
            inx = np.argsort(a.flatten())[::-1]
            nn = int(np.ceil(fraction * a.size))
            if selectfullrank:
                if nn > ngoodrank and ngoodrank > 1:
                    v1 = a[inx[ngoodrank - 1]]
                    v2 = a[inx[ngoodrank]]
                    if verbose:
                        print(f"  good rank threshold: {v1:e} -> {v2:e}")
                        print("  good rank threshold: C {:e} -> {:e}".format(v1**m, v2**m))
            if selectfullrank:
                if verbose:
                    print(" selectfullrank: reducing %d arrays to %d" % (nn, ngoodrank))
                nn = min(nn, ngoodrank)
            if verbose:
                print("final round %d columns: input %d arrays, reduced to %d" % (k, data.size, nn))

            if a.size > 0:
                drawAvalues(a[inx[0:nn]], fig=200 + k)
                plt.title(r"Case $%s$: %d columns, selection, sorted" % (adata.latexstr(), k))

            ABhelper.selectArrays(afile, outfile, inx[0:nn], cache=cache)
            # outfilerest = 'selecteddynamicrest-' + adata0.idstr() + '.oa'
            # ABhelper.selectArrays(afile, outfilerest, inx[nn:], cache=cache)
            cmdlogfid.write("# Select best %d arrays, write to file %s\n\n" % (nn, outfile))

        kn = k + 1
        adatan = oalib.arraydata_t(adata, kn)
        nextfile = "selecteddynamicresult-" + adatan.idstr() + ".oa"
        logfile = nextfile.replace(".oa", "-log.txt")
        dynamicExtendFile(
            outfile, nextfile, kfinal=kfinal, Afinal=Afinal, verbose=1, cmdlog=cmdlogfid, logfile=logfile, cache=cache
        )
        dt = parseProcessingTime(logfile)
        totaltime += dt
        print(f"processing: final step {dt:.0f} [s], {float(dt) / 3600:.1f} [h]")

    afile = nextfile
    adatan = oalib.arraydata_t(adata, kfinal)
    anafiles = analyseFile(afile, method="full", verbose=1, cache=cache)

    anafile = anafiles[1]
    data = loadAnalysisFile(anafile, ncolshint=4)
    A = data[:, 1]
    B = data[:, 2]
    rnk = data[:, 0]
    Acc = 1e-15 ** (1.0 / m)
    idx = A > Acc
    a = A[idx]
    b = B[idx]
    k = adatan.ncols
    m = 1 + k + k * (k - 1) / 2

    nmaxrnk = (rnk == m).nonzero()[0].size
    # Final selection ##

    gidxmask = A >= Afinal

    if selectfullrank:
        print("selectfullrank: go! ")
        gidxmask = (A >= Afinal) * (rnk == m)
    bidxmask = gidxmask is False
    gidx = gidxmask.nonzero()[0]
    bidxmask.nonzero()[0]

    outfile = "finalselection-" + adatan.idstr() + "-A%.2f" % Afinal + ".oa"

    if gidx.size > 300:
        afmode = oalib.ABINARY
    else:
        afmode = oalib.ATEXT

    selectArrays(afile, outfile, gidx, verbose=1, cache=cache, afmode=afmode)

    print("generateABcase: make VIF selection (non-generate)")
    outfileVIF = "finalselection-VIF-" + adatan.idstr() + "-A%.2f" % Afinal + ".oa"

    Bc = B.copy()
    Bc[Bc == 0] = np.inf
    nmaxvif = min(80, B.nonzero()[0].size)
    selectionArraysScore(afile, outfileVIF, Bc, nmax=nmaxvif, order="ascending", cache=cache)

    if a.size > 0:
        print("%d/%d arrays above threshold %.3f. max A value %.5f" % (gidx.size, A.size, Afinal, a.max()))
    else:
        print("%d/%d arrays above threshold %.3f. max A value %.5f" % (gidx.size, A.size, Afinal, 0))

    rr = dict(
        {
            "adata": adata,
            "maxDeff": safemax(a),
            "minB": safemin(b, np.inf),
            "totaltime": totaltime,
            "nselected": gidx.size,
            "narrays": A.size,
            "outfile": outfile,
            "datadir": datadir,
        }
    )
    rr["outfileVIF"] = outfileVIF
    rr["abscatterplot"] = None
    rr["calcmode"] = selectfullrank

    figfile = "finalDAscatter.png"
    if pfig is not None:

        plotABfigure(a, b, figid=100 + 1, verbose=1, fontsize=13)
        plotABboundaries(Acc=Acc)
        ss = adata.latexstr()
        ss = ss.replace("\\OA", r"\mathrm{OA}")
        plt.title(r"Scatterplot for $%s$" % ss, fontsize=18)
        plt.savefig(figfile)

        if verbose:
            print("generating scatterplot: %s" % figfile)
        if gidx.size > 20 and adata.strength < 4 and nmaxrnk > 0:
            rr["abscatterplot"] = figfile

    rr["datadir"] = datadir
    rr["fullcalc"] = selectfullrank == 0 and Afinal == 0
    sys.stdout.flush()
    if Afinal == 0:
        rr["mode"] = "full"
    else:
        rr["mode"] = "partial"

    rr["specialcase"] = 0

    rr["kstart"] = kstart
    rr["Afinal"] = Afinal
    rr["kfinal"] = kfinal
    return rr
    # Estimate total number of arrays


# %%


def tickfontsize(fontsize=14, ax=None):
    """Set the font size for ticks on the x- and y-axis"""
    if ax is None:
        ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    plt.draw()


# %%


def generateStatistics(afile, nbest=10, verbose=1, doprojef=0):
    # FIXME: incomplete function
    #    data=...
    #    gma=...
    data = []
    ll = oalib.readarrayfile(afile)
    data = np.zeros((len(ll), 5))
    for ii, al in enumerate(ll):
        data[ii, 1] = al.Defficiency()
        data[ii, 2] = al.VIFefficiency()
        data[ii, 3] = al.Eefficiency()
        alxf = oalib.array2xf(al)
        data[ii, 0] = alxf.rank()
    gma = np.zeros((len(ll), ll[0].n_columns + 1))
    for ii, al in enumerate(ll):
        gma[:, :] = al.GWLP()
    #        data[:,:]=...

    return generateStatisticsX(afile, data, gma, nbest=nbest, verbose=1)


def generateStatisticsX(afile, data, gma, nbest=10, verbose=1):
    """Convert data and gma to statistics dict"""
    if data.shape[0] <= nbest:
        inx = np.argsort(2 * data[:, 0].flatten() + data[:, 1].flatten())[::-1]
        gmadata = loadGMAdata(afile, cache=1, verbose=0)
        if inx.size == 0:
            # no arrays...
            gsub = gmadata
            ddx = data
        else:
            gsub = gmadata[inx, :]
            ddx = data[inx, :]
        fullstatistics = dict({"mode": "all", "aidx": inx, "data": ddx, "gma": gma, "gwlpdata": gsub})
    else:
        # sort on 2*rank + Defficiency, this is equal to the Defficiency, but
        # numerically more stable
        inx = np.argsort(2 * data[:, 0].flatten() + data[:, 1].flatten())[::-1]
        nn = np.min((nbest, data.shape[0]))
        inx = inx[0:nn]
        gmadata = loadGMAdata(afile, cache=1, verbose=0)
        gsub = gmadata[inx, :]
        # pdb.set_trace()
        # ww=dict({'fulldata': gsub})
        fullstatistics = dict({"mode": "selection", "aidx": inx, "data": data[inx, :], "gma": gma, "gwlpdata": gsub})
        # FIXME: select best arrays
        # rr['fullstatistics']=None
    return fullstatistics


def selectBestArrays(datadir, afile0, outputfile, nmax=None, cache=1, verbose=0):
    """Select arrays from file base on Deffciency score"""
    afile = os.path.join(datadir, afile0)
    data = analyseFileCached(afile, method="rankvalues", verbose=1, cache=cache)
    Deff = data[:, 1]
    rnk = data[:, 0]
    inx = np.argsort(2 * rnk.flatten() + Deff.flatten())[::-1]
    # gc.collect(); print('collect 1 hack'); return dict()

    if nmax is None:
        nn = data.shape[0]
    else:
        nn = np.min((nmax, data.shape[0]))
    if verbose:
        print("selectBestArrays nmax %s, %s: nn %d" % (str(nmax), data.shape, nn))
    # pdb.set_trace()

    # gc.collect(); print('collect 2 hack'); return dict()
    axx = researchOA.makeArraySelection(afile, inx[0:nn], outputfile, verbose=0, dogz=None, cache=cache, hack=0)
    # return dict(); # HACK
    return dict({"aax": axx, "dstr": "%d %s" % (axx[2], axx[1]), "filename": axx[0]})


def analyseABcase(rr, verbose=1, cache=1, subdir="tpages"):
    """Helper function"""
    datadir = rr["datadir"]
    outfile = rr["outfile"]
    rr["arrayfile"] = rr["outfile"]
    ad = rr["adata"]
    t = ad.strength

    if outfile is not None:
        gma = analyseGMA(os.path.join(datadir, outfile), cache=cache, verbose=0)
    else:
        gma = dict({"n": None, "indices": None, "gwp": None})
        # gma=dict({'idx': ii, 'seq': bgma, 'array': None})

    rr["gma"] = gma
    gmadata = rr["gma"]
    # htmldir=rr['htmldir']

    # analyse and select good D-efficiency and VIF values
    afile = os.path.join(rr["datadir"], outfile)
    if ("narrays" in list(rr.keys())) is False:
        rr["narrays"] = nArrayFile(afile)
        rr["nselected"] = nArrayFile(afile)
    if ("totaltime" in list(rr.keys())) is False:
        rr["totaltime"] = None

    rr["ABselection"] = None
    # adatan=oalib.arraydata_t(ad, kfinal)
    anafiles = analyseFile(afile, method="full", verbose=1, cache=cache)
    anafile = anafiles[1]
    data = loadAnalysisFile(anafile, verbose=0)
    Deff = data[:, 1]
    B = data[:, 2]
    Eeff = data[:, 3]
    rnk = data[:, 0]
    k = ad.ncols
    m = 1 + k + k * (k - 1) / 2
    nmaxrnk = (rnk == m).nonzero()[0].size
    ngoodrank = nmaxrnk
    rr["rmaxrnk"] = nmaxrnk
    rr["ngoodrank"] = nmaxrnk
    #    print('x %s' % str(rnk))
    #    print('Deff %s' % str(Deff))
    if verbose:
        print("analyseABcase: k %d: m %d, rmaxrnk %d, ngoodrank %d" % (k, m, nmaxrnk, ngoodrank))

    rr["fullstatistics"] = generateStatisticsX(afile, data, gma)

    if 1:
        if verbose >= 3:
            print("Eff:")
            print(Eeff)
        maxE = ABhelper.safemax(Eeff)
        rr["maxE"] = maxE

        maxDeff = ABhelper.safemax(Deff)
        rr["maxDeff"] = maxDeff

        b = B.copy()
        b[b == 0] = np.inf
        minB = ABhelper.safemin(b, np.inf)
        rr["minB"] = minB

        adelta = 0.00001
        bdelta = 0.00001
        aidx = Deff > (maxDeff - adelta)
        aidx = aidx.nonzero()[0]
        bidx = b <= (minB + bdelta)
        bidx = bidx.nonzero()[0]

        selDfile0 = "selection-goodDeff.oa"
        selBfile0 = "selection-goodB.oa"
        selAfile = os.path.join(datadir, selDfile0)
        selBfile = os.path.join(datadir, selBfile0)
        axx = researchOA.makeArraySelection(afile, aidx, selAfile, verbose=0, dogz=None, cache=cache)
        bxx = researchOA.makeArraySelection(afile, bidx, selBfile, verbose=0, dogz=None, cache=cache)
        if verbose >= 2:
            print("analyseABcase: writing %s " % selBfile)

        if axx[3]:
            selDfile0 += ".gz"
            print("  selAfile: %s adding .gz" % selAfile)
        if bxx[3]:
            selBfile0 += ".gz"

        if verbose >= 2:
            print("  writing AB files %s (%d), %s (%d)" % (selAfile, aidx.size, selBfile, bidx.size))
        rr["ABselection"] = dict({"nA": aidx.size, "nB": bidx.size, "Afile": selDfile0, "Bfile": selBfile0})
        if verbose:
            print(
                "analyseABcase: selected %d good D-efficiency arrays, %d good A-efficiency arrays"
                % (aidx.size, bidx.size)
            )

    # select arrays with good GWLP
    if gmadata["indices"] is not None:
        bgma = np.around(gmadata["gwp"], decimals=6)

        # write array file
        gfile0 = "gmaarray-t%d-%s.oa" % (t, ad.idstr())
        gfile = os.path.join(datadir, gfile0)
        # select the gma arrays

        afile = os.path.join(rr["datadir"], outfile)
        xx = makeArraySelection(afile, gmadata["indices"], gfile, verbose=0, dogz=None, cache=cache)
        pxx = xx[1]
        na = xx[2]

        gstrnolink = ",".join(["%.4f" % v for v in bgma])
        gstr = ",".join(["%.4f" % v for v in bgma])
        # gstr += '  ' + e.small(e.a('%d %s' % (na, pxx), href=subdir+'/' +
        # htmlsubdir0 + '/' + gfile0) )
        gstr += "  " + e.small(e.a("%d %s" % (na, pxx), href=gfile0))
        rr["gmafile"] = gfile0
    else:
        gstr = "-"
        gstrnolink = "-"
        rr["gmafile"] = None
    rr["gmastr"] = gstr
    rr["gmastrnolink"] = gstrnolink

    if ("nselected" in list(rr.keys())) is False:
        rr["nselected"] = -1

    if verbose:
        print(
            "case %s: %d/%d narrays, max D-efficiency %.3f, maxrank %d"
            % (ad.idstr(), rr["nselected"], rr["narrays"], rr["maxDeff"], rr["rmaxrnk"])
        )
    return rr


def copyOAfile(source, targetdir, target0, convert=None, zipfile=None, verbose=1, cache=1):
    """Copy an OA file, depending on arguments convert or compress the file

    Args:
        source (str): source file
        target (str): target output file
        convert (None, int or str): if None copy file. If a number, then convert to either text or binary. If a string, then 'T', 'B', or 'D'
        zipfile (None or bool): of True, compress the file with gzip
    Returns:
        target0final (str): output filename (without .gz extension)
    """
    if convert is None:
        targetfile = os.path.join(targetdir, target0)
        target0final = target0
        if not checkFiles(targetfile, cache):
            if verbose:
                print(f"copyfile {source} -> {targetfile}")
            shutil.copyfile(source, targetfile)
    else:
        na = oapackage.nArrayFile(source)
        if not (isinstance(convert, bytes) or isinstance(convert, str)):
            if na < convert:
                convert = "T"
            else:
                convert = "B"
        if zipfile is None:
            zipfile = convert == "B"
        else:
            if zipfile is not False:
                # print('na %d, zipfile %d' % (na, zipfile))
                zipfile = na >= zipfile
        if not (
            convert == "TEXT"
            or convert == "BINARY"
            or convert == "B"
            or convert == "T"
            or convert == "D"
            or convert == "Z"
            or convert == "DIFF"
        ):
            raise NameError("copyOAfile: convert: should be T, B or D")
        if verbose >= 3:
            print(f"target0: {target0}, zipfile {zipfile}")
        if zipfile:
            if verbose:
                print("copyOAfile: converting to format %s" % convert)
            if target0.endswith(".gz"):
                target0final = target0
                target0 = target0final[:-3]
            else:
                target0final = target0 + ".gz"
            targetfilefinal = os.path.join(targetdir, target0final)
            targetfile = os.path.join(targetdir, target0)
        else:
            target0final = target0
            if target0final.endswith(".gz"):
                print("error: target file ends with .gz")
                raise
            targetfile = os.path.join(targetdir, target0)
            targetfilefinal = os.path.join(targetdir, target0)
        if verbose:
            print(f"copyOAfile: target0 {target0} -> {target0final} ")
        if verbose >= 2:
            print("copyOAfile: converting %s to %s (%d arrays, zip %d)" % (source, targetfile, na, zipfile))
        if checkFiles(targetfilefinal, cache):
            print("  copyOAfile: target file already exist")
        else:
            if verbose >= 2:
                print(f"cmd: oaconvert -v 0 -f {convert} {source} {targetfile}")
            os.system(f"oaconvert -v 0 -f {convert} {source} {targetfile}")
            if zipfile:
                os.system("gzip -f %s" % targetfile)
    return target0final


def convertOAfile(f, verbose=1, dorun=0, maxsize=1000 * 1000, tmpdir=None):
    """Convert OA file to binary compressed form"""
    if not os.path.exists(f):
        if verbose:
            print("file %s does not exist" % f)
        return 0
    if not f.endswith(".oa.gz") and not f.endswith(".oa"):
        if verbose:
            print("file %s does not end with" % f)
        return 0
    sz = os.path.getsize(f)
    if sz > maxsize:
        if verbose:
            print("file %s is too large (%d/%d kb)" % (f, sz / 1024, maxsize / 1024))
        return 0

    af = oalib.arrayfile_t(f, 0)
    isc = af.iscompressed and not af.mode == oalib.AERROR
    af.closefile()
    del af

    if tmpdir is None:
        tmpfile = "tmp.oa"
    else:
        tmpfile = os.path.join(tmpdir, "tmp.oa")
    tmpfilegz = tmpfile + ".gz"

    if isc:
        if verbose >= 2:
            print("file %s is already compressed" % f)
        return 0

    if f.endswith(".oa.gz"):
        gzx = 1
        fplain = f.replace(".oa.gz", ".oa")
        if os.path.exists(fplain):
            if verbose:
                print("file %s is already there?!?" % fplain)
            return 0
    else:
        gzx = 0
        fplain = f

    if gzx:
        cmd = "gzip -d %s" % f
    else:
        cmd = 'echo ""'
    # cmd2='oajoin -i %s -f B -o tmp' % fplain
    cmd2 = f"oaconvert -f B {fplain} {tmpfile}"
    cmd3 = "gzip %s" % tmpfile
    cmd4 = f"rm -f {fplain}; mv {tmpfilegz} {f}"
    cmd5 = "oainfo %s" % f
    if verbose >= 2:
        print("Commands:")
        print(cmd)
        print(cmd2)
        print(cmd3)
        print(cmd4)
        print(cmd5)
    if dorun:
        os.system(cmd)
        os.system(cmd2)
        os.system(cmd3)
        os.system(cmd4)
        os.system(cmd5)
    return 1


#    oadir = fullfile(basedir, sprintf('design-%s', fullstr));
#    return oadir
def defaultCase(N, t, kmax, s, casedir=None):
    s += [s[-1]] * (kmax - len(s))
    si = oalib.intVector(s)
    ad = oalib.arraydata_t(si, N, t, kmax)
    return dict({"case": ad, "casedir": caseDefaultDir(ad)})


def readNumbersFile(nfile, verbose=1):
    fid = open(nfile)
    # d = fid.readline()
    nums = []
    # print('d: %s' % d)
    while 1:
        d = fid.readline()
        # if verbose>=2:
        #    print(d)
        if verbose >= 2:
            print("line %s" % d)
        if not d[0:2] == "k ":
            break
        if verbose >= 2:
            print("line " + d)
        kn = scanf.sscanf(d, "k %d: %d %d")
        nums.append(kn)
    fid.close()
    return nums


def caseReadNumbers(case, basedatadir, verbose=1):
    ad = case["case"]
    numbersfile = "numbers-%s-t%d.txt" % (ad.idstr(), ad.strength)
    xf = os.path.join(basedatadir, case["casedir"], numbersfile)
    fid = open(xf)
    d = fid.readline()
    nums = []
    while 1:
        d = fid.readline()
        if verbose >= 2:
            print("line %s" % d)
        if not d[0:2] == "k ":
            break
        if verbose >= 2:
            print("line " + d)
        kn = scanf.sscanf(d, "k %d: %d arrays")
        nums.append(kn[1])
    fid.close()
    return nums


# 3

# %%


def evenoddClusterGetkmax(N, strength):
    if strength == 3:
        k = int(np.ceil(1 + N / 2))
    else:
        k = int(N - 1)
    if k > 30:
        k = 30
    if N == 48 and strength == 3:
        k = 18
    if N == 56 and strength == 3:
        # special case to make equal to N=64 case
        k = 24
    if N == 64 and strength == 3:
        # special case to reduce the total number of files
        k = 24
    return k


# %%


def evenoddCases(N, strength=3, lastlevel=3):
    if N == 32:
        splitdata = dict()
        splitdata["N"] = N
        splitdata["strength"] = strength
        splitdata["kinitial"] = 6
        splitdata["klevel"] = 7
        splitdata["klevel2"] = 10
        splitdata["klevelcheck"] = 10
        splitdata[0] = dict({"n": 2, "k": splitdata["kinitial"]})
        splitdata[1] = dict({"n": 2, "k": splitdata["klevel"]})
        iisel = list(range(0, splitdata[0]["n"]))
        jjsel = list(range(0, splitdata[1]["n"]))
    if N == 40:
        splitdata = dict()
        splitdata["N"] = N
        splitdata["strength"] = strength
        splitdata["kinitial"] = 6
        splitdata["klevel"] = 7
        splitdata["klevel2"] = 8
        splitdata["klevelcheck"] = 10
        splitdata[0] = dict({"n": 2, "k": splitdata["kinitial"]})
        splitdata[1] = dict({"n": 4, "k": splitdata["klevel"]})
        iisel = list(range(0, splitdata[0]["n"]))
        jjsel = list(range(0, splitdata[1]["n"]))
    if N == 48:
        splitdata = dict()
        splitdata["N"] = N
        splitdata["strength"] = strength
        splitdata["kinitial"] = 6
        splitdata["klevel"] = 7
        splitdata["klevelcheck"] = 12
        splitdata["klevel2"] = 9
        splitdata[0] = dict({"n": 2, "k": splitdata["kinitial"]})
        splitdata[1] = dict({"n": 4, "k": splitdata["klevel"]})
        iisel = list(range(0, splitdata[0]["n"]))
        jjsel = list(range(0, splitdata[1]["n"]))

    if N == 56:
        lastlevel = 2
        splitdata = dict()
        splitdata["N"] = N
        splitdata["strength"] = strength
        splitdata["kinitial"] = 7
        splitdata["klevel"] = 9
        splitdata["klevel2"] = 9
        splitdata["klevelcheck"] = 12
        splitdata[0] = dict({"n": 10, "k": splitdata["kinitial"]})
        splitdata[1] = dict({"n": 12, "k": splitdata["klevel"]})
        iisel = list(range(0, splitdata[0]["n"]))
        jjsel = list(range(0, splitdata[1]["n"]))
    if N == 64:
        splitdata = dict()
        splitdata["N"] = N
        splitdata["strength"] = strength
        splitdata["klevelcheck"] = 14
        splitdata["kinitial"] = 7
        splitdata["klevel"] = 9
        splitdata["klevel2"] = 10  # special case
        splitdata[0] = dict({"n": 780, "k": splitdata["kinitial"]})
        splitdata[1] = dict({"n": 248, "k": splitdata["klevel"]})
        iisel = list(range(0, splitdata[0]["n"]))
        jjsel = list(range(0, splitdata[1]["n"]))

    splitdata["kmax"] = evenoddClusterGetkmax(N, strength)
    if N == 40 and strength == 3:
        splitdata["kmax"] = 14
    #    if N==64 and strength==3:
    #        splitdata['kmax']=24
    if N == 48:
        splitdata[2] = dict({"n": 8, "k": splitdata["klevel2"]})
    else:
        splitdata[2] = dict({"n": 18, "k": splitdata["klevel2"]})
        splitdata[3] = dict({"n": 36, "k": splitdata["klevel2"] + 1})

    splitdata["lastlevel"] = lastlevel
    splitdata["levels"] = [5] + [splitdata[i]["k"] for i in range(lastlevel)] + [splitdata["kmax"]]

    return (splitdata, iisel, jjsel)


# %%


def jobStatus(alljobs, verbose=1):
    """Print status for a list of jobs"""
    print("job status: gathered %d jobs" % len(alljobs))

    jobs = [j for j in alljobs if not j.complete()]
    gjobs = [j for j in jobs if j.canrun()]
    if verbose >= 2:
        for i, j in enumerate(jobs):
            print("job %d: %s" % (i, j))
    if verbose:
        print("  %d/%d to run" % (len(gjobs), len(jobs)))

    return gjobs, jobs


def processingtimeFile(lvls):
    return "processing-%s.txt" % splitTag(lvls)


def mdFile(lvls):
    mdfileX = "md5check-%s.txt" % splitTag(lvls)
    return mdfileX


def numbersFile(lvls, tag="numbers"):
    return f"{tag}-{splitTag(lvls)}.txt"


def writeprocessingTime(pfile, dt):
    try:
        f = open(pfile, "w")
        f.write("%f\n" % dt)
        f.close()
    except:
        raise
        pass
    return True


# def readProcessingFile(fname):
#    try:
#        f=open(fname, 'rt')
#        x=f.readline().strip()
#        T=float(x)
#        f.close()
#        return T
#    except:
#        return None


def readprocessingTime(pfile):
    """Read processing time from a file"""
    try:
        f = open(pfile)
        s = f.readline().strip()  # = f.read()
        f.close()
        T = float(s)
    except:
        T = None
        pass
    return T


# %%


def evenoddAnalyseRun(outputdir, adata, splitdata, verbose=1, iimax=None):
    ptimes = []
    badidx = []
    nbx = []

    if iimax is None:
        iimax = splitdata[0]["n"]
    for ii in range(0, iimax):
        edir = splitdir([ii])
        nfile = os.path.join(outputdir, edir, "timelog-split-%d.txt" % ii)
        try:
            fid = open(nfile)
        except:
            badidx.append(ii)
            continue
        jj = 0
        ptsub = []
        for ln in fid:
            if ln.startswith("#time total"):
                dt = scanf.sscanf(ln, "#time total: %f [s]")
                jj = jj + 1
                if verbose >= 2:
                    print("time %.2f " % dt)
                ptimes.append((ii, jj, dt[0]))
                ptsub.append((ii, jj, dt[0]))
        fid.close()
        ttsub = [x[2] for x in ptsub]
        if verbose:
            if len(ttsub) > 0:
                print(
                    "evenoddAnalyseRun: processing time block %d: %.1f [h], mean: %.1f +- %.1f [s]"
                    % (ii, np.sum(ttsub) / 3600.0, np.mean(ttsub), np.std(ttsub))
                )

    tt = [x[2] for x in ptimes]

    for kk in range(5, adata.ncols):
        adata0 = oapackage.arraydata_t(adata, kk)
        outfile = os.path.join(outputdir, "results-j5evenodd-pareto-%s.oa" % adata0.idstr())
        if verbose:
            na = oapackage.nArrays(outfile)
            print("  column %d: %d Pareto designs" % (kk, na))

    if len(tt) > 0:
        print(
            "evenoddAnalyseRun: processing time: %.1f [h], mean: %.1f +- %.1f [s]"
            % (np.sum(tt) / 3600, np.mean(tt), np.std(tt))
        )
        print(
            " gathered %d/%d files (fraction %.3f)"
            % (len(tt), splitdata[0]["n"] * splitdata[1]["n"], float(len(tt)) / (splitdata[0]["n"] * splitdata[1]["n"]))
        )

        ncores = 100
        th = np.sum(tt) / 3600
        tEstimate = (1.0 / 24) * (1.0 / float(ncores)) * th * ((splitdata[0]["n"] * splitdata[1]["n"]) / len(tt))
        print("evenoddAnalyseRun: time estimated: cluster %.1f [d] (with %d cores)" % (tEstimate, ncores))

    return (ptimes, badidx, nbx)


def doSplitFile(lvls, splitdata, adata, verbose=1, outputdir=None, createdirs=False, cache=True):
    splittag = ".".join([str(x) for x in lvls])
    level = len(lvls)
    rfile = splitname(lvls)
    if verbose:
        print(f"split level {splittag}: file {rfile} (extend and split)")
    ebase = rfile.replace(".oa", "-extend")
    edir = splitdir(lvls)
    mkdirc(edir)

    rfile.replace(".oa", "-sp%d" % level)
    nn = splitdata[level]["n"]
    splitname(lvls + [nn - 1])
    splitout = rfile.replace(".oa", "-sp%d" % (level))

    klevel = splitdata[level]["k"]
    adatax = oapackage.arraydata_t(adata, klevel)

    if splitdata["levels"][level] == splitdata["levels"][level + 1]:
        # special case: in the previous step no extension was performed
        edirup = splitdir(lvls[:-1])
        afile = os.path.join(edirup, splitname(lvls))
    else:
        afile = os.path.join(edir, "{}-{}".format(ebase, adatax.idstr() + ".oa"))

    if verbose:
        print(f"split level {splittag}: file {afile}")
    splitcmd = "oasplit -v 1 -f Z --nb 1 -i %s -n %d -o %s; " % (afile, nn, os.path.join(edir, splitout))
    # cmd += splitcmd
    done = False
    checkfile = splitname(lvls + [splitdata[level]["n"] - 1])
    # pdb.set_trace()

    if outputdir is not None:
        # pdb.set_trace()
        if createdirs:
            for kk in range(splitdata[level]["n"]):
                x = splitdir(lvls + [kk])
                if verbose:
                    print("doSplitFile: make dir %s" % x)
                _ = oapackage.mkdirc(join(outputdir, x))

        if oahelper.checkArrayFile(os.path.join(outputdir, os.path.join(edir, checkfile)), cache=cache):
            done = True
    return splitcmd, afile, splitout, done


# %%


def checkLevel(lvls, splitdata, adata, outputdir, verbose=1):
    """Check whether all data at a certain level has been generated"""
    edir = splitdir(lvls)
    tag = splitTag(lvls)
    # check numbers file
    level = len(lvls)
    if verbose:
        print("checkLevel: %s" % tag)
    # check Pareto files
    kmin = splitdata["levels"][level] + 1
    kmid = splitdata["levels"][level + 1]
    kmax = splitdata["kmax"]

    fulldata = 1
    if verbose:
        print("checkLevel: %s: column %d to %d and %d to %d" % (tag, kmin, kmid, kmid + 1, kmax))
    for ii, k in enumerate(range(kmin, kmid + 1)):
        adatax = oapackage.arraydata_t(adata, k)
        vv = adatax.idstr()

        afile0 = splitBase(lvls) + "-extend-%s.oa" % vv
        afile = os.path.join(outputdir, edir, afile0)
        if not oahelper.checkFilesOA(afile):
            if verbose:
                print(f"checkLevel {tag}: no result file ({afile0})")
            fulldata = 0
            break
    for ii, k in enumerate(range(kmid + 1, kmax + 1)):
        adatax = oapackage.arraydata_t(adata, k)
        vv = adatax.idstr()
        afile0 = splitBase(lvls) + "-pareto-%s.oa" % vv
        afile = os.path.join(outputdir, edir, afile0)
        if not oahelper.checkFilesOA(afile):
            if verbose:
                print(f"checkLevel {tag}: no result file ({afile0})")
            fulldata = 0
            break

    if fulldata == 0:
        if verbose:
            print("checkLevel %s: no result data present (k %d)" % (tag, k))
        return 0

    # check numbers
    nfile = numbersFile(lvls)  # 'numbers-%s.txt' % tag
    numbersfile = os.path.join(outputdir, edir, nfile)
    if not oahelper.checkFiles(numbersfile):
        if verbose:
            print(f"checkLevel {tag}: no numbers file ({numbersfile})")
        return 0
    return 1


# %%


def touchFile(fname, txt=""):
    with open(fname, "w") as f:
        f.write(txt)
    return


def pythonBinary():
    if platform.system() == "Linux":
        return "python3"
    return "python"


def isVSCcluster():
    if "VSC_SCRATCH" in os.environ.keys():
        return True
    else:
        return False


def check_pipe_cmd(tag="oaclustergather", errorcode=1):
    cmd = "\n"
    cmd += "RES=${PIPESTATUS[0]}" + os.linesep
    cmd += "if [ $RES -eq 0 ]" + os.linesep
    cmd += "then" + os.linesep
    cmd += '  echo "   successfully executed %s"' % tag + os.linesep
    cmd += "else" + os.linesep
    cmd += '  echo "   problem with %s: errorcode $?"' % tag + os.linesep
    cmd += "  exit %d" % errorcode + os.linesep
    cmd += "fi" + os.linesep
    return cmd


def gatherLockfile(outputdir, lvls, tag=None):
    edir0 = splitdir(lvls)
    if tag is None:
        tag = splitTag(lvls)
    lockfile0 = "lockfile-%s.txt" % tag
    lockfile = join(outputdir, edir0, lockfile0)
    return lockfile


def gatherResults(
    lvls,
    outputdir,
    splitdata,
    adata,
    dozip=True,
    verbose=2,
    nparetodiff=1,
    paretomethod=1,
    gatherpareto=True,
    gatherJstats=True,
    legacy=False,
    ncores=1,
    queue="q72h",
    quiet=False,
):
    """Gather results from multiple jobs

    lvls : list
        level for gather
    legacy : bool
        parameter passed to gatherFilesList
    dozip (bool):
        zip results into a single zip file
    """
    joblist = []

    tag = splitTag(lvls)
    edir0 = splitdir(lvls)
    level = len(lvls)
    kmin = splitdata["levels"][level + 1] + 1
    kmax = splitdata["kmax"]
    prfile = os.path.join(outputdir, edir0, processingtimeFile(lvls))
    nfile = os.path.join(outputdir, edir0, numbersFile(lvls))
    nfilej = os.path.join(outputdir, edir0, numbersFile(lvls, tag="jstats"))
    zipfile0 = "splitted-%s.zip" % tag
    zipfile = join(outputdir, edir0, zipfile0)
    lockfile = gatherLockfile(outputdir, lvls)
    plog = os.path.join(outputdir, edir0, "gather-log-%s.txt" % tag)
    jlog = os.path.join(outputdir, edir0, "gather-log-jstats-%s.txt" % tag)

    if os.path.exists(lockfile):
        print("gatherResults: lockfile %s exists" % lockfile)
        return []
    pythonbinary = pythonBinary()
    if len(lvls) > 0:
        ii = lvls[0]
    else:
        ii = -1
    if len(lvls) > 1:
        jj = lvls[1]
    else:
        jj = -1
    if len(lvls) > 2:
        kk = lvls[2]
    else:
        kk = -1
    nsplitstr = " ".join(["--nsplit%d %d" % (x, splitdata[x]["n"]) for x in range(4) if x in splitdata])
    cmd = f'echo "Calculating pareto optimal arrays"; cd {join(outputdir, edir0)};\n'
    cmd += ("touch %s" % lockfile) + os.linesep

    lvloptions = " --split0 %d --split1 %d --split2 %d %s --kmax %d --kmin %d" % (ii, jj, kk, nsplitstr, kmax, kmin)
    if gatherpareto:
        cmd += (
            "oaclustergather -c %s -b %s -f BINARY --paretomethod %d --nparetodiff %d --numbersfile %s %s --cleanrun 1 -o %s | tee %s;"
            % (
                join(outputdir, "oaconfig.txt"),
                join(outputdir),
                paretomethod,
                nparetodiff,
                nfile,
                lvloptions,
                "pareto",
                plog,
            )
        )
        cmd += os.linesep
    if gatherJstats:
        cmd += os.linesep + "# j-statistics " + os.linesep
        cmd += "oaclustergather -c {} -b {} --method 1 --numbersfile {} {} --cleanrun 1 -o {} | tee {};".format(
            join(outputdir, "oaconfig.txt"),
            join(outputdir),
            nfilej,
            lvloptions,
            "jstats",
            jlog,
        )
        cmd += os.linesep
    if isVSCcluster():
        cmd += check_pipe_cmd(tag="oaclustergather", errorcode=1)

    cmd += os.linesep + 'echo "paretomethod %d" > %s\n' % (paretomethod, join(outputdir, edir0, "gatheroptions.txt"))
    cmd += gzipOA(join(outputdir, edir0))
    cmd += os.linesep
    cmd += "cd %s; tail -n 3 sp%d-split-*/sp0*extend.txt > %s;" % (
        join(outputdir, edir0),
        level,
        os.path.join("timelog-split-%s.txt" % tag),
    )
    cmd += "cd %s;" % outputdir
    cmd += os.linesep

    # FIXME: hard-coded
    codedir = os.path.join(os.path.expanduser("~"), "projects/oapackage")
    if not os.path.exists(codedir):
        raise Exception("update hard-coded location!")
    gfile = os.path.join(codedir, "pythondevelop", "gather_processing_times.py")
    cmd += "\n%s %s -b %s -n %d %s" % (
        pythonbinary,
        gfile,
        outputdir,
        splitdata[level]["n"],
        " ".join([str(x) for x in lvls]),
    )
    cmd += os.linesep

    if dozip:
        zipstr = "".join(["sp%d-split*" % i for i in range(level + 1)]) + ".oa*"
        # zipstr = '-'.join(['sp%d-split-%d'  % (i, lvls[i]) for i in range(level) ]) + '-sp%d-split-*.oa*'  % (level+1)

        cmdz = "\n# compress subdirs" + os.linesep
        cmdz += ("PFILE=%s" % prfile) + os.linesep
        cmdz += "if [ -f $PFILE ];" + os.linesep
        cmdz += "then" + os.linesep
        if quiet:
            qstr = " -q "
        else:
            qstr = ""
        cmdz += (
            'echo "Compressing directory %s"; cd %s; zip -r -m %s -T %s sp%d-split-*  %s'
            % (os.path.join(outputdir, edir0), os.path.join(outputdir, edir0), qstr, zipfile, level, zipstr)
            + os.linesep
        )
        cmdz += os.linesep
        cmdz += "else" + os.linesep
        cmdz += ('   echo "Could not find file %s"' % prfile) + os.linesep

        cmdz += "fi" + os.linesep

        cmd += cmdz
    else:
        cmd += os.linesep + "# no compressing of subdirs (dozip=False)" + os.linesep

    cmd += ("rm -f %s" % lockfile) + os.linesep
    # runcommand(cmd, dryrun=dryrun, idstr='job-compress-%d (depends on job-%d-x and job-%d-statistics)' % (ii,ii, ii), verbose=verbose, logfile=allcmdlogfile)

    for kk in range(splitdata[level]["n"]):
        lvlsn = lvls + [kk]
        edirX = splitdir(lvlsn)

        lfile = os.path.join(os.path.join(outputdir, edirX), "lock-extend.txt")

        if os.path.exists(lfile):
            print(f"gatherResults {tag}: lock file for {kk} exists, skipping")
            return joblist

    checkfilesstart = gatherFilesList(lvls, outputdir, splitdata, adata, verbose=0, legacy=legacy)
    if dozip:
        checkfiles = [nfile, prfile, zipfile]
    else:
        checkfiles = [nfile, prfile]
    j = job(
        cmd=cmd,
        jobtype="gather %s" % tag,
        queue=queue,
        shorttag="G%s" % tag,
        ncores=ncores,
        checkfiles=checkfiles,
        checkfilesstart=checkfilesstart,
    )

    if not j.complete() or True:
        joblist.append(j)
    return joblist


# %%
def gatherFilesList(lvls, outputdir, splitdata, adata, verbose=1, paretofiles=True, legacy=False):
    """Return list of files needed to gather results for given stage"""
    tag = splitTag(lvls)
    splitdir(lvls)
    level = len(lvls)
    kmin = splitdata["levels"][level + 1] + 1
    kmax = splitdata["levels"][level + 2]
    # nfile = os.path.join(outputdir, edir0, numbersFile(lvls))
    # flist=[nfile]
    flist = []
    kmin = splitdata["levels"][level + 1] + 1
    kmid = splitdata["levels"][level + 2]
    kmax = splitdata["kmax"]

    lastlevel = splitdata["lastlevel"] == (level + 1)

    if verbose:
        print("gatherFiles: %s: kmin %d, kmid %d, kmax %d, lastlevel %d" % (tag, kmin, kmid, kmax, lastlevel))

    for kk in range(splitdata[level]["n"]):
        lvlsn = lvls + [kk]
        # tagX = splitTag(lvlsn)
        edir = splitdir(lvlsn)
        for ii, k in enumerate(range(kmin, kmid + 1)):
            adatax = oapackage.arraydata_t(adata, k)
            vv = adatax.idstr()

            afile0 = splitBase(lvls + [kk]) + "-extend-%s.oa" % vv
            afile = os.path.join(outputdir, edir, afile0)
            flist += [afile]

        for ii, k in enumerate(range(kmid + 1, kmax + 1)):
            adatax = oapackage.arraydata_t(adata, k)
            vv = adatax.idstr()

            afile0 = splitBase(lvls + [kk]) + "-pareto-%s.oa" % vv
            afile = os.path.join(outputdir, edir, afile0)
            if not legacy:
                flist += [afile]

        if kmid < kmax:
            if not legacy:
                pfile = os.path.join(outputdir, edir, processingtimeFile(lvlsn))
                flist += [pfile]
        if lastlevel:
            pfile = os.path.join(outputdir, edir, mdFile(lvlsn))

            flist += [pfile]

    return flist


def evenoddAnalysePartialRun(outputdir, adata, splitdata, verbose=1, iimax=None):
    adataxx = oapackage.arraydata_t(adata, splitdata["klevelcheck"])

    ntotal = 0
    nfiles = 0
    totaltime = 0
    if iimax is None:
        iimax = splitdata[0]["n"]
    for ii in range(0, iimax):
        if ii % 80 == 0:
            print("evenoddAnalysePartialRun: parsing %d/%d" % (ii, splitdata[0]["n"]))
        edir0 = splitdir([ii])
        pfile = researchOA.processingtimeFile([ii])
        T = readprocessingTime(join(outputdir, edir0, pfile))

        nfile = researchOA.numbersFile([ii])

        if T is not None:
            # this level is already complete
            totaltime += T
            nfiles += splitdata[1]["n"]
            aa = researchOA.loadnumbersfile(join(outputdir, edir0, nfile))
            ntotal += aa[aa[:, 0] == splitdata["klevelcheck"], :][0, 1]
            if verbose >= 2:
                print("evenoddAnalysePartialRun: file %d/%d already done" % (ii, splitdata[0]["n"]))

            continue

        for jj in range(0, min(splitdata[1]["n"], 100000)):

            rfile = splitname([ii, jj])
            ebase = rfile.replace(".oa", "-extend")
            edir = splitdir([ii, jj])

            afile = os.path.join(edir, "{}-{}".format(ebase, adataxx.idstr() + ".oa"))

            cmdlogfile = os.path.join(edir, rfile.replace(".oa", "-extend.txt"))
            dt0 = oapackage.parseProcessingTime(cmdlogfile, verbose=0)
            if dt0 >= 0:
                totaltime += dt0

            v = oapackage.nArrays(afile)
            if v >= 0:
                if verbose >= 2:
                    print("file %d-%d: %d (time %.1f [s])" % (ii, jj, v, dt0))
                elif verbose >= 1 and v > 10:
                    print("file %d-%d: %d (time %.1f [s])" % (ii, jj, v, dt0))

                ntotal += v
                nfiles += 1
    nnn = splitdata[0]["n"] * splitdata[1]["n"]
    if nfiles == 0:
        sfac = 1
    else:
        sfac = float(splitdata[0]["n"] * splitdata[1]["n"]) / float(nfiles)

    ncores = 100
    print(
        "evenoddAnalysePartialRun: %d/%d files extended, # arrays at col %d: %d (estimate %.2e)"
        % (nfiles, nnn, splitdata["klevelcheck"], ntotal, 4 * ntotal * sfac)
    )
    print("  time {:.1f} [s], estimate {:.1f} [h] ".format(totaltime, float(totaltime) * sfac / 3600.0))
    print(
        "  time estimated: cluster %.1f [d] (with %d cores)"
        % (float(totaltime) * sfac / (ncores * 24.0 * 3600.0), ncores)
    )

    xfac = splitdata[0]["n"] * splitdata[1]["n"]
    edir = splitdir([0, 3])
    if os.path.exists(edir):
        cfactor = arrayfilecompression(edir, verbose=0)
        print(
            "expected disk usage: %.3f GB, %.1f bytes per array"
            % (xfac * float(get_folder_size(edir)) / 1024.0**3, cfactor)
        )

    return (ntotal, nfiles)


def casesOA(allcases=1):
    """Return all cases with results"""
    cases = []

    cases.append(defaultCase(4, 2, 4, [2, 2, 2, 2]))
    cases.append(defaultCase(8, 2, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(8, 2, 6, [4, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(9, 2, 5, [3, 3, 3, 3, 3]))
        cases.append(defaultCase(12, 2, 4, [6, 2, 2, 2]))
        cases.append(defaultCase(12, 2, 6, [3, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(12, 2, 12, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(16, 2, 10, [8, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(16, 2, 6, [4, 4, 4, 4, 4, 4]))
        cases.append(defaultCase(16, 2, 8, [4, 4, 4, 4, 2, 2, 2, 2]))
        cases.append(defaultCase(16, 2, 10, [4, 4, 4, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(16, 2, 12, [4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(16, 2, 14, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(16, 2, 16, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(18, 2, 7, [6, 3, 3, 3, 3, 3, 3]))
        cases.append(defaultCase(18, 2, 8, [3, 3, 3, 3, 3, 3, 3, 3]))
        cases.append(defaultCase(18, 2, 9, [2, 3, 3, 3, 3, 3, 3, 3, 3]))
        cases.append(defaultCase(20, 2, 4, [10, 2, 2, 2]))
        cases.append(defaultCase(20, 2, 10, [5, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(20, 2, 20, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(20, 2, 10, [5, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(24, 2, 14, [12, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-24-122a-t2"

    cases.append(defaultCase(24, 2, 20, [6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-24-62a"

    cases.append(defaultCase(24, 2, 22, [4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-24-432a"

    if allcases:
        cases.append(defaultCase(24, 2, 14, [6, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(24, 2, 23, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-24-42a-t2"
        cases.append(defaultCase(24, 2, 20, [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-24-32a-t2"

    cases.append(defaultCase(24, 2, 24, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-24-2a-t2"

    if allcases:
        cases.append(defaultCase(25, 2, 6, [5, 5, 5, 5, 5, 5]))
        cases.append(defaultCase(27, 2, 12, [9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        cases.append(defaultCase(27, 2, 14, [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        cases[-1]["casedir"] = "design-27-3a"

    cases.append(defaultCase(28, 2, 14, [7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-28-72a-t2"

    cases.append(
        defaultCase(28, 2, 27, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    )
    cases[-1]["casedir"] = "design-28-2a-t2"

    cases.append(defaultCase(30, 2, 3, [5, 3, 2]))
    cases.append(defaultCase(32, 2, 20, [2]))
    cases.append(defaultCase(32, 2, 10, [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]))
    cases[-1]["casedir"] = "design-32-4a-t2"

    if allcases:
        cases.append(defaultCase(32, 2, 22, [4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-4444442a-t2"

        cases.append(defaultCase(32, 2, 18, [4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-44444442a-t2"

        cases.append(defaultCase(32, 2, 18, [4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-444444442a-t2"

        cases.append(defaultCase(32, 2, 18, [4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-4444444442a-t2"

        cases.append(defaultCase(32, 2, 14, [8, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-84444442a-t2"

        cases.append(defaultCase(32, 2, 16, [8, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-844444442a-t2"

        cases.append(defaultCase(32, 2, 20, [8, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-8444442a-t2"

        cases.append(defaultCase(32, 2, 20, [8, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-32-844442a-t2"

        cases.append(defaultCase(N=32, t=2, kmax=18, s=[4, 4, 4, 2]))
    cases.append(defaultCase(N=32, t=2, kmax=18, s=[4, 4, 2]))
    cases.append(defaultCase(N=32, t=2, kmax=18, s=[4, 2]))

    if allcases:
        cases.append(defaultCase(36, 2, 6, [6, 6, 6, 6, 6, 6]))
        cases[-1]["casedir"] = "design-36-6a-t2"

        cases.append(defaultCase(36, 2, 17, [6, 6, 6, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        cases[-1]["casedir"] = "design-36-6663a-t2"

        cases.append(defaultCase(36, 2, 16, [6, 6, 6, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-36-66632a-t2"

        cases.append(defaultCase(36, 2, 20, [6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-36-6662a-t2"

        cases.append(defaultCase(42, 2, 3, [7, 3, 2]))
        cases.append(defaultCase(70, 2, 3, [7, 5, 2]))
        cases[-1]["casedir"] = "design-70-752a-t2"

    cases.append(defaultCase(36, 2, 10, [2, 2]))
    cases.append(defaultCase(N=36, t=2, kmax=8, s=[3, 2]))
    cases.append(defaultCase(N=36, t=2, kmax=8, s=[3, 3, 2]))

    cases.append(defaultCase(N=45, t=2, kmax=8, s=[5, 3]))
    cases.append(defaultCase(N=108, t=3, kmax=8, s=[4, 3]))
    cases.append(defaultCase(N=60, t=2, kmax=8, s=[5, 3, 2]))
    cases.append(defaultCase(N=40, t=2, kmax=8, s=[5, 4, 2]))
    cases.append(defaultCase(N=40, t=2, kmax=8, s=[4, 2]))
    cases.append(defaultCase(N=40, t=2, kmax=8, s=[5, 2]))

    ############## strength 3 #################################
    cases.append(defaultCase(8, 3, 5, [2, 2, 2, 2, 2]))
    cases.append(defaultCase(16, 3, 5, [4, 2, 2, 2, 2]))
    cases.append(defaultCase(16, 3, 9, [2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(24, 3, 5, [6, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(24, 3, 6, [3, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(24, 3, 13, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(27, 3, 5, [3, 3, 3, 3, 3]))
        cases.append(defaultCase(32, 3, 5, [8, 2, 2, 2, 2]))
        cases.append(defaultCase(32, 3, 7, [4, 4, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(32, 3, 9, [4, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(32, 3, 17, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(36, 3, 4, [3, 3, 2, 2]))
    cases.append(defaultCase(N=36, t=2, kmax=18, s=[3]))
    cases.append(defaultCase(40, 3, 6, [10, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-40-102a"

    if allcases:
        cases.append(defaultCase(40, 3, 8, [5, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(40, 3, 21, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(48, 3, 5, [12, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-48-122a-t3"

    cases.append(defaultCase(48, 3, 5, [6, 4, 2, 2, 2]))
    cases.append(defaultCase(48, 3, 9, [6, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(48, 3, 17, [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-48-32a"

    cases.append(defaultCase(48, 3, 7, [4, 3, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(48, 3, 14, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-48-42a-t3"

    cases.append(defaultCase(48, 3, 25, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-48-2a-t3"

    if allcases:
        cases.append(defaultCase(54, 3, 4, [6, 3, 3, 3]))
        cases.append(defaultCase(54, 3, 7, [2, 3, 3, 3, 3, 3, 3]))
        cases.append(defaultCase(54, 3, 6, [3, 3, 3, 3, 3, 3]))

    cases.append(defaultCase(56, 3, 5, [14, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-56-142a-t3"

    cases.append(defaultCase(56, 3, 7, [7, 2, 2, 2, 2, 2, 2]))
    cases.append(
        defaultCase(
            56, 3, 32, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        )
    )
    cases[-1]["casedir"] = "design-56-2a-t3"

    cases.append(defaultCase(60, 3, 4, [5, 3, 2, 2]))
    cases.append(defaultCase(64, 3, 5, [16, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-64-162a-t3"

    cases.append(defaultCase(64, 3, 4, [8, 4, 2, 2]))
    cases.append(defaultCase(64, 3, 9, [8, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-64-82"

    cases.append(defaultCase(64, 3, 7, [4, 4, 4, 4, 4, 4, 4]))
    cases.append(defaultCase(64, 3, 8, [4, 4, 4, 4, 4, 4, 2, 2]))
    if allcases:
        cases.append(defaultCase(64, 3, 9, [4, 4, 4, 4, 4, 2, 2, 2, 2]))
        cases.append(defaultCase(64, 3, 11, [4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(64, 3, 13, [4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(64, 3, 18, [4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-64-442a-t3"

    cases.append(defaultCase(64, 3, 17, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-64-42a-t3"

    cases.append(
        defaultCase(
            64,
            3,
            35,
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
        )
    )
    cases[-1]["casedir"] = "design-64-2a-t3"

    if allcases:
        cases.append(defaultCase(72, 3, 4, [4, 3, 3, 2]))
        cases.append(defaultCase(72, 3, 15, [3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-72-332a-t3"

    cases.append(defaultCase(72, 3, 9, [9, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-72-92a-t3"

    if allcases:
        cases.append(defaultCase(72, 3, 5, [6, 6, 2, 2, 2]))
        cases.append(defaultCase(72, 3, 9, [6, 3, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(80, 3, 10, [5, 4, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-80-542a-t3"

    cases.append(defaultCase(80, 3, 10, [2, 2, 2, 2, 2, 2, 2, 2]))

    cases.append(defaultCase(81, 3, 6, [9, 3, 3, 3, 3, 3]))
    cases[-1]["casedir"] = "design-81-93a-t3"

    cases.append(defaultCase(81, 3, 12, [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
    cases[-1]["casedir"] = "design-81-3a-t3"

    if allcases:
        cases.append(defaultCase(84, 3, 4, [7, 3, 2, 2]))
        cases.append(defaultCase(88, 3, 14, [2]))
        cases.append(defaultCase(90, 3, 4, [5, 3, 3, 2]))
        cases.append(defaultCase(96, 3, 14, [2]))

        cases.append(defaultCase(96, 3, 14, [3, 2]))

        cases.append(defaultCase(100, 3, 4, [5, 5, 2, 2]))
        cases.append(defaultCase(104, 3, 14, [2]))

        cases.append(defaultCase(108, 3, 4, [3, 3, 2, 2]))
        cases[-1]["casedir"] = "design-108-3322-t3"

        cases.append(defaultCase(108, 3, 5, [3, 3, 3, 2, 2]))
        cases[-1]["casedir"] = "design-108-33322-t3"

        cases.append(defaultCase(108, 3, 4, [6, 3, 3, 2]))
        cases[-1]["casedir"] = "design-108-6332-t3"

        cases.append(defaultCase(112, 3, 4, [7, 4, 2, 2]))
        cases[-1]["casedir"] = "design-112-7422-t3"

        cases.append(defaultCase(120, 3, 4, [5, 4, 3, 2]))
        cases[-1]["casedir"] = "design-120-5432-t3"

        cases.append(defaultCase(120, 3, 4, [6, 5, 2, 2]))
        cases[-1]["casedir"] = "design-120-6522-t3"

        cases.append(defaultCase(150, 3, 4, [5, 5, 3, 2]))
        cases[-1]["casedir"] = "design-150-5532a-t3"

        cases.append(defaultCase(176, 3, 14, [2]))

    cases.append(defaultCase(16, 4, 6, [2, 2, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(32, 4, 6, [4, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(32, 4, 7, [2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(48, 4, 5, [6, 2, 2, 2, 2]))
        cases.append(defaultCase(48, 4, 7, [3, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(48, 4, 6, [2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(64, 4, 5, [8, 2, 2, 2, 2]))
        cases.append(defaultCase(64, 4, 5, [4, 4, 2, 2, 2]))
        cases.append(defaultCase(64, 4, 8, [4, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(64, 4, 9, [2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(72, 4, 5, [3, 3, 2, 2, 2]))
    cases.append(defaultCase(80, 4, 6, [5, 2, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(80, 4, 6, [2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(81, 4, 5, [3, 3, 3, 3, 3]))
        cases.append(defaultCase(96, 4, 6, [6, 4, 2, 2, 2, 2]))
        cases.append(defaultCase(96, 4, 7, [6, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(96, 4, 6, [4, 3, 2, 2, 2, 2]))
        cases.append(defaultCase(96, 4, 8, [4, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(96, 4, 9, [3, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(96, 4, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(108, 4, 5, [3, 3, 3, 2, 2]))
        cases.append(defaultCase(112, 4, 6, [7, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(112, 4, 6, [2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(120, 4, 5, [5, 3, 2, 2, 2]))
    cases.append(defaultCase(128, 4, 6, [16, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-128-162a-t4"

    cases.append(defaultCase(128, 4, 12, [8, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-128-82a-t4"

    cases.append(defaultCase(128, 4, 12, [8, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-128-842a"

    cases.append(defaultCase(128, 4, 8, [4, 4, 4, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(128, 4, 9, [4, 4, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(128, 4, 16, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-128-42a"

    cases.append(defaultCase(128, 4, 18, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-128-2a"

    cases.append(defaultCase(144, 4, 11, [9, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-144-92a"

    cases.append(defaultCase(144, 4, 11, [6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-144-662a"

    if allcases:
        cases.append(defaultCase(144, 4, 6, [6, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(144, 4, 6, [6, 3, 2, 2, 2, 2]))
        cases.append(defaultCase(144, 4, 5, [4, 3, 3, 2, 2]))
        cases.append(defaultCase(144, 4, 10, [3, 3, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-144-332a-t4"

        cases.append(defaultCase(144, 4, 8, [3, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-144-32a"

        cases.append(defaultCase(144, 4, 12, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-144-2a"

        cases.append(defaultCase(160, 4, 7, [10, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-160-102a-t4"

        cases.append(defaultCase(160, 4, 8, [5, 4, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-160-542a-t4"

        cases.append(defaultCase(160, 4, 11, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-160-42a-t4"

        cases.append(defaultCase(160, 4, 10, [5, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-160-52a-t4"

    cases.append(defaultCase(160, 4, 11, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-160-2a-t4"

    if allcases:
        cases.append(defaultCase(162, 4, 6, [3, 3, 3, 3, 3, 3]))
        cases.append(defaultCase(162, 4, 6, [3, 3, 3, 3, 3, 2]))
        cases.append(defaultCase(162, 4, 5, [3, 3, 3, 3, 2]))
        cases.append(defaultCase(168, 4, 5, [7, 3, 2, 2, 2]))
        cases.append(defaultCase(176, 4, 14, [2]))
        cases.append(defaultCase(180, 4, 5, [5, 3, 3, 2, 2]))
        cases.append(defaultCase(200, 4, 5, [5, 5, 2, 2, 2]))
        cases.append(defaultCase(216, 4, 5, [4, 3, 3, 3, 2]))
        cases.append(defaultCase(216, 4, 5, [6, 3, 3, 2, 2]))
        cases.append(defaultCase(252, 4, 5, [7, 3, 3, 2, 2]))
        cases.append(defaultCase(256, 4, 6, [4, 4, 4, 4, 4, 4]))
        cases[-1]["casedir"] = "design-256-4a-t4"

        cases.append(defaultCase(216, 4, 6, [3, 3, 3, 2, 2, 2]))
        cases[-1]["casedir"] = "design-216-3332a-t4"

        cases.append(defaultCase(256, 4, 14, [4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-256-444442a-t4"

    # strength 5 and higher
    # cases.append(defaultCase(16, 4, 5, [ 2,2,2,2,2 ]))
    cases.append(defaultCase(32, 5, 7, [2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(64, 5, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(64, 5, 7, [4, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(96, 5, 7, [2, 2, 2, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(96, 5, 8, [3, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(128, 5, 10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(128, 5, 8, [4, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(128, 5, 7, [4, 4, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(128, 5, 7, [8, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(144, 5, 6, [3, 3, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(160, 5, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(160, 5, 8, [5, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 5, 8, [4, 3, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 5, 9, [2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 5, 7, [4, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 5, 8, [6, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 5, 9, [3, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 5, 7, [6, 4, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(216, 5, 6, [3, 3, 3, 2, 2, 2]))
        cases.append(defaultCase(224, 5, 8, [7, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-224-72a-t5"

    if allcases:
        cases.append(defaultCase(224, 5, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(240, 5, 6, [5, 3, 2, 2, 2, 2]))
        cases.append(defaultCase(256, 5, 9, [4, 4, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(256, 5, 7, [4, 4, 4, 2, 2, 2, 2]))
        cases.append(defaultCase(288, 5, 6, [4, 3, 3, 2, 2, 2]))
        cases.append(defaultCase(360, 5, 6, [5, 3, 3, 2, 2, 2]))
        cases.append(defaultCase(384, 5, 15, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-384-42a-t5"

    cases.append(defaultCase(64, 6, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(128, 6, 9, [2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(128, 6, 8, [4, 2, 2, 2, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(192, 6, 9, [3, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(192, 6, 10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(256, 6, 8, [4, 4, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(256, 6, 9, [4, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases.append(defaultCase(256, 6, 10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-256-2a-t6"

    cases.append(defaultCase(288, 6, 7, [3, 3, 2, 2, 2, 2, 2]))
    if allcases:
        cases.append(defaultCase(320, 6, 8, [2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(384, 6, 12, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-384-2a-t6"

    cases.append(defaultCase(384, 6, 12, [3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-384-32a-t6"

    if allcases:
        cases.append(defaultCase(384, 6, 10, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-384-42a-t6"

    cases.append(defaultCase(384, 6, 11, [4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-384-432a-t6"

    if allcases:
        cases.append(defaultCase(128, 7, 9, [2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(256, 7, 10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-256-2a-t7"

        cases.append(defaultCase(256, 7, 10, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-256-42a-t7"

        cases.append(defaultCase(384, 7, 10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-384-2a-t7"

        cases.append(defaultCase(384, 7, 10, [3, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-384-32a-t7"

        cases.append(defaultCase(384, 7, 9, [6, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(512, 7, 13, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-512-2a-t7"

        cases.append(defaultCase(512, 7, 11, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-512-42a-t7"

        cases.append(defaultCase(512, 7, 10, [4, 4, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-512-442a-t7"

        cases.append(defaultCase(256, 8, 10, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases.append(defaultCase(512, 8, 12, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
        cases[-1]["casedir"] = "design-512-2a-t8"

    cases.append(defaultCase(512, 8, 10, [4, 2, 2, 2, 2, 2, 2, 2, 2, 2]))
    cases[-1]["casedir"] = "design-512-42a-t8"
    return cases
