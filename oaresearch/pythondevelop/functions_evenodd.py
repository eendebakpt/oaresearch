"""
Created on Fri Nov 16 11:24:39 2012

@author: eendebakpt
"""


""" Load necessary packages """


import datetime
import math
import os
import time
from os.path import join

import numpy as np
import oalib
import oapackage
from oapackage import oahelper
from oaresearch.pythondevelop.researchOA import *
from oaresearch.pythondevelop.researchOA import (
    analyseFile,
    checkFiles,
    checkFilesOA,
    evenoddCases,
    formatProccessingTime,
    gatherLockfile,
    gatherResults,
    loadAnalysisFile,
    nArrayFile,
    parseProcessingTime,
    readNumbersFile,
    runcommand,
    splitdir,
    splitname,
    splitTag,
)

import researchOA
from ABhelper import *


def smallEOcase(key: str, adata, basedir, htmldir: str, latexdir: str, subdir=None, verbose=1, cache=1):
    """Calculate one of the small cases"""
    t = adata.strength
    if subdir is None:
        subdir = "evenodd-" + adata.idstr() + "-t%d" % adata.strength
        casedir = os.path.join(basedir, subdir)
        if not os.path.exists(casedir):
            os.mkdir(casedir)
        adata.writeConfigFile(os.path.join(casedir, "oaconfig.txt"))

        cmd = "cd %s; oaextendsingle -f B -l 2  | tee log0.txt " % casedir
        if verbose >= 2:
            print(cmd)
        rfile = os.path.join(casedir, "result-%s.oa" % adata.idstr())
        if checkFiles(rfile, cache=cache):
            if verbose >= 2:
                print("file already exists")
        else:
            os.system(cmd)
    else:
        casedir = os.path.join(basedir, subdir)

    rr = dict({"adata": adata, "casedir": casedir})
    rr["totaltime"] = parseProcessingTime(os.path.join(casedir, "log0.txt"))
    rr["subdir"] = subdir

    with open(os.path.join(casedir, "code.txt"), "wt") as fid:
        fid.write("generated with smallEOcase\n")
        fid.write("paretomethod: 1\n")

    adata = rr["adata"]
    nb = adata.ncols + 2
    maxrnk = np.zeros((adata.ncols + 1,))
    eotable = np.zeros((t, nb))
    eotable[t - 1, 0] = 1
    par = dict()
    maxrnk[t] = 1 + t + t * (t - 1) / 2  # the root is the only possible array
    for ii in range(t + 1, adata.ncols + 1):
        adata0 = oalib.arraydata_t(adata, ii)
        afile = os.path.join(casedir, "result-%s.oa" % adata0.idstr())
        print(afile)
        eotype = geteotype(afile, cache=cache, verbose=0)
        ww = np.histogram(eotype, range(0, nb + 1))
        if verbose >= 4:
            print("%s" % str(ww))
        if verbose >= 1:
            print("  kk %d: %d arrays %s" % (ii, nArrayFile(afile), ww[0]))
        eotable = np.vstack((eotable, ww[0]))

        pfile = os.path.join(casedir, "result-pareto-%s.oa" % adata0.idstr())
        cmd = "oapareto %s --paretomethod 1 -f T -v 1 -o %s" % (afile, pfile)
        if checkFiles(pfile, cache=cache):
            if verbose >= 2:
                print("file %s already exists" % pfile)
        else:
            os.system(cmd)
        nn = oalib.nArrays(pfile)
        adata0 = oalib.arraydata_t(rr["adata"], int(ii))
        pfile = os.path.join(casedir, "result-pareto-%s.oa" % adata0.idstr())
        sols = oalib.readarrayfile(pfile)
        xx = oalib.parsePareto(sols, 0)
        if verbose:
            print("pareto: %d classes, %d elements" % (xx.number(), xx.numberindices()))
        par[ii] = dict({"narrays": nn, "nclasses": xx.number()})

        for jj, al in enumerate(sols):
            r = oalib.array2xf(al).rank()
            if maxrnk[ii] < r:
                maxrnk[ii] = r

    kmax = eotable.sum(axis=0).nonzero()[0][-1]
    ncolsmax = eotable.sum(axis=1).nonzero()[0][-1] + 1
    rr["kmax"] = kmax
    rr["ncolsmax"] = ncolsmax
    rr["pareto"] = par
    rr["maxrank"] = maxrnk
    print("  kmax %d" % kmax)
    rr["eotable"] = eotable
    rr["narrays"] = eotable.sum(axis=1)
    rr["eoseries"] = eotable[:, 1:].sum(axis=1)
    rr["evenseries"] = eotable[:, 0:].sum(axis=1)
    rr["subpage"] = None

    filename = os.path.join(latexdir, "eotable-%s.txt" % adata.fullidstr(series=1))

    writeEOtable(adata, eotable, filename, t=t, maxrnk=maxrnk, kmax=kmax, verbose=2)
    print("  written to file %s" % filename)
    return rr


def evenoddSubpageSpecial(htmldir, r, verbose=1, pagecache=1):
    """Create html page for even-odd cases calculated on the cluster"""
    casedir = r["casedir"]
    adata = r["adata"]
    xstr = series2htmlstr(adata, case=0)
    xstrplain = series2htmlstr(adata, html=0, case=0)
    subfile = "evenodd-%s-t%d.html" % (adata.idstr(), adata.strength)
    strength = adata.strength

    page = markup.page()
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
                "description": "Even-Odd arrays",
            }
        ),
        footer="<!-- End of page -->",
    )

    page.h1("Series %s " % xstr)
    oap = e.a("Orthogonal Array", href="../software.html")
    pstr = "This page contains information about even-odd and foldover arrays. "
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

    tstart = int(1 + 2 * np.floor(float(strength) / 2))

    page.table()
    page.tr(style="font-weight: bold;")

    if 0:
        eotable = r["eotable"]
        kmax = eotable.sum(axis=0).nonzero()[0][-1]

        ss = ["Class", "Foldover"] + ["E%d" % x for x in range(tstart, kmax + 3, 2)]
        page.td(ss[0], style="padding-right:30px;")
        page.td(ss[1:], style="padding-right:8px;")
    #        page.tr.close()
    else:
        ss = ["Class", "Number of arrays"]  # (of J5 type)
        page.td(ss[0:], style="padding-right:8px;")

    page.td("Pareto arrays (pareto classes)", style="padding-right:8px;")
    page.td("Max rank /  Max D-efficiency", style="padding-right:8px;")

    page.tr.close()

    kinitial = r["kinitial"]
    nn = adata.ncols
    #  nn=eotable.sum(axis=1).nonzero()[0][-1]+1
    # nn=adata.ncols

    nnp = r["nnp"]
    for ii in range(adata.strength, nn + 1):
        adata0 = oalib.arraydata_t(adata, ii)
        xstr0 = series2htmlstr(adata0, case=1)

        pfile0 = "results-j5evenodd-pareto-%s.oa" % adata0.idstr()
        pfile = os.path.join(casedir, pfile0)
        if ii in nnp.keys():
            na = nnp[ii]["narrays"]
            nt = nnp[ii].get("ntotal", -1)
            # npar=nnp[ii]['np']
            if nt > -1:
                nfo = nnp[ii]["nfoldover"]
                neo = nt - nnp[ii]["nfoldover"]
                ss = ["total %d (foldover %d, even-odd %d)" % (nt, nfo, neo)]
            else:
                ss = ["j5-type %d" % (na)]
        else:
            ss = ["?"]
            na = -1

        pnstr = "Series: " + e.a(xstr, href=subfile) + ", " + e.a("main page", href="../index-evenodd.html")
        # generate subpage
        if ii >= 5 and na != 0:
            dstr = "Pareto optimal designs of J5 type." % ()
            targetdir = os.path.join(htmldir, "evenoddpages")
            subfilebase = "paretoclassx-%s" % adata0.fullidstr()
            subfilecase = parseParetoClass(
                adata0, pfile, dstr=dstr, targetdir=targetdir, subfilebase=subfilebase, prevnextstr=pnstr, verbose=1
            )
            xstrlnk = e.a(xstr0, href=subfilecase)
        else:
            xstrlnk = xstr0

        page.tr(style="")
        page.td(xstrlnk, style="padding-right:30px;")

        page.td(ss, style="padding-right:8px;")

        parstr = "%d (%d)" % (nnp[ii]["nparetoarrays"], nnp[ii]["nparetoclasses"])
        page.td(parstr, style="padding-right:8px;")

        m = 1 + ii + ii * (ii - 1) / 2
        if nnp[ii]["maxrnk"] > -1:
            if nnp[ii]["maxrnk"] == m:
                parstr = "<b>%d</b> / %.3f" % (nnp[ii]["maxrnk"], nnp[ii]["maxDeff"])
            else:
                parstr = "%d / %.3f" % (nnp[ii]["maxrnk"], nnp[ii]["maxDeff"])
        else:
            parstr = "-"

        page.td(parstr, style="padding-right:8px;")
        page.tr.close()

        if na == 0 and nt == -1:
            break

    page.table.close()

    page.br(clear="both")
    tstr = "The max rank is the maximum rank of the second order interaction matrix of the arrays."
    page.p(tstr)

    tstr = e.a("main page", href="../index-evenodd.html")
    page.p(tstr)

    if "totaltime" in r:
        page.br(clear="both")
        tstr = formatProccessingTime(r["totaltime"])
        page.p("Processing time: &nbsp; %s" % tstr)
    #  page.p( e.small('kmax %d, adata.ncols %d, nn %d' % (kmax, adata.ncols, nn)) )

    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)
    page.br(clear="both")
    page.p("Page generated on %s." % dstr)

    r["subpage0"] = subfile
    fid = open(os.path.join(htmldir, "evenoddpages", subfile), "wt")
    fid.write(str(page))
    fid.close()
    return subfile, adata


def evenoddSubpage(htmldir, r, verbose=1, pagecache=1):
    casedir = r["casedir"]
    adata = r["adata"]
    t = adata.strength
    xstr = series2htmlstr(adata, case=0)
    xstrplain = series2htmlstr(adata, html=0, case=0)
    subfile = "evenodd-%s-t%d.html" % (adata.idstr(), adata.strength)
    strength = adata.strength

    page = markup.page()
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
                "description": "Even-Odd arrays",
            }
        ),
        footer="<!-- End of page -->",
    )

    page.h1("Series %s " % xstr)
    oap = e.a("Orthogonal Array", href="../software.html")
    pstr = "This page contains information about even-odd and foldover arrays. "
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

    tstart = int(1 + 2 * np.floor(float(t) / 2))

    page.table()
    page.tr(style="font-weight: bold;")

    eotable = r["eotable"]
    kmax = eotable.sum(axis=0).nonzero()[0][-1]

    ss = ["Class", "Foldover"] + ["E%d" % x for x in range(tstart, kmax + 3, 2)]
    page.td(ss[0], style="padding-right:30px;")
    page.td(ss[1:], style="padding-right:8px;")
    page.td("Total", style="padding-right:8px; font-weight: bold; ")
    page.td("Max rank", style="padding-right:8px;")
    page.tr.close()

    nn = r["ncolsmax"] + 1
    #  nn=eotable.sum(axis=1).nonzero()[0][-1]+1
    # nn=adata.ncols

    maxrnk = r["maxrank"]
    for ii in range(adata.strength, nn):
        adata0 = oalib.arraydata_t(adata, ii)
        xstr0 = series2htmlstr(adata0, case=1)

        print("evenoddSubpage: columns %d, kmax: %s, eotable.shape %s" % (ii, kmax, eotable.shape))
        vals = np.hstack((eotable[ii - 1, 0:1], eotable[ii - 1, tstart : (kmax + 3) : 2]))
        if verbose >= 2:
            print(vals)
        ss = []  # [' %d  ' % ( ii )]
        ss += ["%d" % x for x in list(vals)]

        nt = sum(vals)
        ss += ["%d" % (sum(vals))]
        m = 1 + ii + ii * (ii - 1) / 2
        if maxrnk[ii] == m and nt > 0:
            tstr = "<b>%d (%d)</b>" % (maxrnk[ii], m)
        else:
            tstr = "%d (%d)" % (maxrnk[ii], m)
        ss += [tstr]

        pnstr = "Series: " + e.a(xstr, href=subfile) + ", " + e.a("main page", href="../index-evenodd.html")
        # generate subpage
        if vals.sum() > 0 and ii > strength:
            dstr = "There are <b>%d</b> non-isomorphic arrays in this class, <b>%d</b> of them are even-odd." % (
                vals.sum(),
                vals[1:].sum(),
            )
            print("evenoddSubpage: vals %d %d" % (vals.sum(), vals[1:].sum()))
            targetdir = os.path.join(htmldir, "evenoddpages")
            pfile = os.path.join(casedir, "result-pareto-%s" % adata0.idstr() + ".oa")
            subfilebase = "paretoclassx-%s" % adata0.fullidstr()
            subfilecase = parseParetoClass(
                adata0, pfile, dstr=dstr, targetdir=targetdir, subfilebase=subfilebase, prevnextstr=pnstr, verbose=1
            )
            xstrlnk = e.a(xstr0, href=subfilecase)
        else:
            xstrlnk = xstr0

        page.tr(style="")
        page.td(xstrlnk, style="padding-right:30px;")

        page.td(ss, style="padding-right:8px;")
        page.tr.close()

    page.table.close()

    page.br(clear="both")
    tstr = "For each number of columns the number of the different types of arrays are printed."
    tstr += " The maximum rank is the maximum of the second order interaction matrix of the arrays. The number between "
    tstr += "parenthesis is the maximum rank that can be achieved."
    page.p(tstr)

    page.br(clear="both")
    tstr = formatProccessingTime(r["totaltime"], keep_seconds=True)
    page.p("Processing time: &nbsp; %s" % tstr)

    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)
    #    page.br(clear='both')
    page.p("Page generated on %s." % dstr)

    page.br(clear="both")
    page.p(e.small("kmax %d, adata.ncols %d, nn %d" % (kmax, adata.ncols, nn)))

    r["subpage0"] = subfile
    fid = open(os.path.join(htmldir, "evenoddpages", subfile), "wt")
    fid.write(str(page))
    fid.close()
    return subfile


def eoseries2html(eotable, r):
    if eotable is None:
        # special case
        adata = r["adata"]
        k = adata.ncols
        ss = ""
        dd = [r["nnp"][ii].get("narrays", -1) for ii in range(adata.strength, k + 1)]
        ddt = [r["nnp"][ii].get("ntotal", -1) for ii in range(adata.strength, k + 1)]
        dd = [-1] * adata.strength + dd
        ddt = [-1] * adata.strength + ddt

        k = min(adata.ncols, (np.array(dd) == 0).nonzero()[0][0] + 1)
        for jj in range(1, k + 1):
            # r['nnp'][jj]
            # ss+=e.span('%d' % eoseries[jj], style='color: blue;')
            if ddt[jj] > -1:
                ss += "(" + e.span("%d" % ddt[jj], style="color: blue;") + ")"
            else:
                if dd[jj] > -1:
                    ss += "(" + e.span("%d" % dd[jj], style="color: orange;") + ")"
                else:
                    ss += "(" + e.span("%d" % 0, style="color: blue;") + ")"

            if jj < 10:
                ss += ", "

        return ss
    ncolsmax = eotable.sum(axis=1).nonzero()[0][-1] + 1
    ncolsmax = min(ncolsmax, eotable.shape[0] - 1)

    eoseries = eotable[:, 1:].sum(axis=1)
    evenseries = eotable[:, 0:].sum(axis=1)
    ii = r["eoseries"].size
    ss = ""
    for jj in range(0, ncolsmax + 1):
        ss += e.span("%d" % eoseries[jj], style="color: blue;")
        ss += "(" + e.span("%d" % evenseries[jj], style="color: red;") + ")"
        if jj < ncolsmax:
            ss += ", "
    return ss


def geteotype(afile, cache=1, verbose=1):
    """Determine even-odd type for all arrays in a file

    Returns:
        eotype (array): array of j-types. for an even array the value is 0.
            For value k we have J_k != 0 and J_l =0 for all odd l < k

    """
    anafile = analyseFile(afile, method="gma", verbose=0, cache=cache)

    if len(anafile) == 0:
        raise Exception("geoeotype: issue with %s" % afile)
    data = loadAnalysisFile(anafile[0], verbose=verbose)

    w = data[:, 1::2]  # contains j1, j3, j5, ...
    eotype = np.zeros(w.shape[0], dtype=int)
    for j in range(0, data.shape[0]):
        if verbose >= 2:
            oapackage.tprint("array %d" % j)
        b = np.nonzero(w[j, 0:] > 1e-10)[0]
        if verbose >= 3:
            print("series %s, b %s" % (w[j, 0:], str(b)))
        if b.size > 0:
            eotype[j] = 1 + 2 * b[0]
        else:
            eotype[j] = 0
    return eotype


def writeEOtable(adata, eotable, filename, t=1, maxrnk=None, kmax=10, verbose=1):
    fid = open(filename, "w")
    tstart = int(1 + 2 * np.floor(float(t) / 2))
    fid.write("%% kmax %d, tstart %d\n" % (kmax, tstart))
    fid.write("%% created by writeEOtable (from functions_evenodd.py)\n")
    v = kmax - tstart + 4
    # fid.write('\\begin{tabular}{%s}\n' % ( 'l'*v ))
    # ss=' & '.join( ['Column', 'Foldover']+['E%d' % x for x in range(tstart,kmax+3,2) ]) + ' \\\\'  + '\n'
    # ss+= '\hline \n'
    # fid.write(ss)
    tmp = eotable.sum(axis=1)
    nonzeroind = np.nonzero(tmp)[0]
    colmax = nonzeroind[-1] + 2
    if verbose:
        print("writeEOtable: %s: colmax %d" % (adata.fullidstr(series=1), colmax))
    Xh = np.array(["Column", "Foldover"] + ["E%d" % x for x in range(tstart, kmax + 3, 2)])
    if maxrnk is not None:
        Xh = np.hstack((Xh, np.array(["Max rank"])))
    X = np.zeros((colmax - t + 1, Xh.shape[0]), dtype=object)
    for ii in range(t, colmax + 1):
        vals = np.hstack((eotable[ii - 1, 0:1], eotable[ii - 1, tstart : (kmax + 3) : 2]))
        if verbose >= 2:
            print("ii: %d, eotableindex %d" % (ii, ii - 1))
            print(vals)
        X[ii - t, 0] = "%d" % ii
        for jj, v in enumerate(vals):
            X[ii - t, jj + 1] = "%d" % v
        ss = " %d   &  " % (ii)
        ss += " & ".join(["%d" % x for x in list(vals)]) + "  \\\\" + "\n"
        # fid.write(ss)
        if maxrnk is not None:
            m = 1 + ii + ii * (ii - 1) / 2
            r = maxrnk[ii]
            if m == r:
                X[ii - t, -1] = "{\\bfseries{%d}}/%d" % (r, m)
            else:
                X[ii - t, -1] = "%d/%d" % (r, m)

    # fid.write('\\end{tabular}\n')
    X = np.vstack((Xh, X))
    ss = researchOA.array2latex(X, header=1, tabchar="l", hlines=[0])
    fid.write(ss)
    fid.close()


# %%


def rerunGather(
    lvls,
    adata,
    splitdata,
    outputdirtmp=None,
    verbose=1,
    cache=1,
    reruntag="rerun",
    gatherpareto=True,
    gatherJstats=True,
    hack=0,
    sublevels=True,
):
    t0 = time.time()
    edir0 = splitdir(lvls)
    tag = splitTag(lvls)
    level = len(lvls)
    splitfile = "splitted-%s.zip" % splitTag(lvls)
    splitfiletag0 = "splitted-%s-%s.txt" % (splitTag(lvls), reruntag)
    splitfiletag = os.path.join(outputdirtmp, edir0, splitfiletag0)

    lockfile = gatherLockfile(outputdirtmp, lvls)

    if os.path.exists(lockfile):
        raise Exception("lockfile %s exists" % lockfile)
    if oahelper.checkFiles(splitfiletag, cache=cache):
        print("rerunGather %s: split directory with tag %s already calculated " % (tag, splitfiletag0))
        return

    if outputdirtmp is None:
        outputdirtmp = tempfile.mkdtemp("oa-%s" % tag)

    if oahelper.checkFiles(os.path.join(outputdirtmp, edir0, splitfile), cache=cache):
        print("rerunGather %s: split directory %s already calculated " % (tag, edir0))
    else:
        if (verbose >= 2) or (verbose >= 1 and level <= 2):
            print("rerunGather %s: split file %s does not exist " % (tag, join(edir0, splitfile)))
        adatax = adata.reduceColumns(splitdata["kmax"])
        edir = splitdir(lvls)
        rfile = splitname(lvls)
        afile = os.path.join(edir, "%s-%s" % (rfile.replace(".oa", "-extend"), adatax.idstr() + ".oa"))
        if oapackage.checkArrayFile(join(outputdirtmp, afile)):
            if (verbose >= 2) or (verbose >= 1 and level <= 2):
                print("   file %s does exist!" % join(afile))
            pass
        else:
            print("error: no final result!!!!")
            exit(1)
        return []

    qqdir = join(outputdirtmp, edir0)

    if sublevels:
        print("rerunGather %s: unzip to %s" % (tag, qqdir))
        cmd = "cd %s; unzip -q -o %s;" % (qqdir, join(outputdirtmp, edir0, splitfile))
        if verbose >= 2:
            print(cmd)
        if verbose:
            pass
            # print(oapackage.__file__)
        r = runcommand(cmd, shell=True)

        if r > 0:
            print("error in rerunGather with zipfile! (errorcode %d)" % r)
            raise Exception("error with zipfile")
            return None
        # r = subprocess.check_output(cmd, shell=True)

        print("rerunGather %s: unzip complete" % (tag,))

        for ij in range(splitdata[level]["n"]):
            sublvls = lvls + [ij]
            if verbose:
                print("rerunGather: sub %s" % (sublvls))

            if hack:
                print("hack: early abort! (no sublevels!!!)")
                break

            rerunGather(
                sublvls,
                adata,
                splitdata,
                outputdirtmp=outputdirtmp,
                reruntag=reruntag,
                gatherpareto=gatherpareto,
                gatherJstats=gatherJstats,
            )
            if len(lvls) <= 1:
                time.sleep(0.002)

            if ij > 40000:
                print("hack: early abort!")
                break

    if splitdata["lastlevel"] == len(lvls):
        print("warning level structure too deep...")
        if lvls[0] == 567 and lvls[1] == 111 and lvls[2] == 16 and adata.N == 64 and adata.strength == 3:
            print("!!! special case, overriding splitdata !!!")
            print("lvls %s" % (lvls,))
            splitdata, iisel, jjsel = evenoddCases(adata.N, strength=adata.strength, lastlevel=4)
        else:
            raise

    if verbose:
        print("rerunGather %s: before gatherResults: %.1f [seconds]" % (tag, time.time() - t0))

    joblistg1 = gatherResults(
        lvls,
        outputdirtmp,
        splitdata,
        adata=adata,
        verbose=1,
        legacy=True,
        quiet=True,
        gatherpareto=gatherpareto,
        gatherJstats=gatherJstats,
    )
    if verbose:
        print("\nrerunGather %s: run gather: %.1f [seconds]" % (tag, time.time() - t0))
        if verbose >= 2:
            print(joblistg1[0].cmd)
            print("---")

    # return None
    if verbose and hack:
        for jx in joblistg1:
            print("rerunGather: level %s: hack: print job command" % tag)
            print(jx.cmd)

        pass
        # print(runcommand.__module__)
        # print(oapackage.__file__)
        # print(cmd)
    for jx in joblistg1:
        jx.execute = True
    jr = [jx.runjob(cache=False) for jx in joblistg1]

    alljobs = joblistg1

    if len(lvls) <= 1:
        time.sleep(0.05)

    good = np.all(jr)
    if good:
        ss = "jobresults: " + str(jr)
        researchOA.touchFile(splitfiletag, "%s\noapackage %s\n%s\n" % (oapackage.timeString(), oapackage.version(), ss))
    else:
        print("rerunGather %s: gather returned error" % (tag,))

    if verbose:
        print("rerunGather %s: done: %.1f [seconds]" % (tag, time.time() - t0))

    return alljobs


# %%


def finalStageProcessing(outputdir, adata, splitdata, lvls, verbose=1, cache=True):
    edir = splitdir(lvls)
    rfile = splitname(lvls)
    tag = splitTag(lvls)

    # make processing and numbers file
    pfile0 = researchOA.processingtimeFile(lvls)
    pfile = os.path.join(outputdir, edir, pfile0)
    if os.path.exists(pfile) and cache:  # and os.path.exists(nfile):
        pass
    else:
        efile = os.path.join(edir, rfile.replace(".oa", "-extend.txt"))
        dtx = oapackage.parseProcessingTime(efile, verbose=0)
        print("generateLevelNext %s: legacy:  read processing time from %s: %.1f" % (tag, efile, dtx))
        researchOA.writeprocessingTime(pfile, dtx)
        # create numbers file?
    return None


def generateLevelNext(
    adata,
    lvls,
    splitdata,
    outputdir,
    cache=True,
    allcmdlogfile=None,
    discardJ5=-1,
    paretomethod=1,
    ncores=None,
    verbose=1,
    splitmode=False,
    priorextend=False,
):
    """Generate jobs scripts for given level of calculation"""
    N = adata.N
    runcomplete = 1
    edir0 = splitdir(lvls[:-1])
    edir = splitdir(lvls)
    rfile = splitname(lvls)
    level = len(lvls)
    tag = splitTag(lvls)
    if ncores is None:
        ncores = 4

    if verbose >= 2:
        print("split tag %s: file %s" % (tag, rfile))
    ebasefull = os.path.join(edir, rfile).replace(".oa", "-extend")
    afilee = os.path.join(edir0, rfile)

    if verbose >= 2:
        print("## generateLevelNext: tag %s" % tag)

    mdfilelegacy = "md5check.txt"
    mdfileX0 = researchOA.mdFile(lvls)

    legacyfile = os.path.join(outputdir, edir, "legacy.txt")
    if checkFiles(legacyfile):
        if verbose >= 2:
            print("generateLevelNext: %s: legacy format (generation complete)" % tag)
        joblist = []
        runcomplete = None
        cmd = None
        return joblist, runcomplete, cmd

    if checkFiles(join(outputdir, edir, mdfilelegacy)):
        # make Pareto files from full extend files
        kmid = splitdata["levels"][level + 1]
        knext = splitdata["levels"][level + 2]

        if verbose:
            print(
                "generateLevelNext: %s: update legacy format: level %d: kmid %d, knext %d: convert extend to pareto"
                % (tag, level, kmid, knext)
            )

        if 1:
            joblist = []
            cc = ""
            for ii, k in enumerate(range(kmid + 1, knext + 1)):
                if verbose >= 2:
                    print("generateLevelNext: %s: update legacy format file: %d " % (tag, k))
                adatax = oapackage.arraydata_t(adata, k)
                vv = adatax.idstr()

                afile0 = splitBase(lvls) + "-extend-%s.oa" % vv
                afile = os.path.join(outputdir, edir, afile0)
                pfile0 = splitBase(lvls) + "-pareto-%s.oa" % vv
                pfile = os.path.join(outputdir, edir, pfile0)
                if os.path.exists(pfile):
                    # print(pfile)
                    continue
                if 0:
                    print("   parse %s to %s" % (afile0, pfile0))
                    researchOA.parseParetoFile(afile, pfile)
                cmd = "oapareto %s --paretomethod %d -o %s" % (afile, paretomethod, pfile)
                cmd += "  | tee -a log-legacy-%s.txt" % tag
                cc += cmd + "\n"
                j = job(
                    cmd,
                    jobtype="convert legacy %s k %d" % (tag, k),
                    shorttag="legacy-%s-%d" % (tag, k),
                    checkfiles=[pfile],
                    checkfilesstart=[afile],
                )
                joblist.append(j)

            if len(joblist) > 0:
                joblist = [
                    job(
                        cc,
                        jobtype="convert legacy %s" % tag,
                        shorttag="legacy-%s" % (tag,),
                        checkfiles=[pfile],
                        checkfilesstart=[afile],
                        ncores=1,
                    )
                ]
                runcomplete = None
                cmd = None
                if verbose:
                    print("generateLevelNext %s: legacy format: return joblist" % tag)
                return joblist, runcomplete, cmd

            # make processing and numbers file
            finalStageProcessing(outputdir, adata, splitdata, lvls, verbose=1, cache=False, paretomethod=paretomethod)

            pfile0 = researchOA.processingtimeFile(lvls)
            nfile0 = researchOA.numbersFile(lvls)
            pfile = os.path.join(outputdir, edir, pfile0)
            nfile = os.path.join(outputdir, edir, nfile0)
            if os.path.exists(pfile):  # and os.path.exists(nfile):
                pass
            else:
                efile = os.path.join(edir, rfile.replace(".oa", "-extend.txt"))
                dtx = oapackage.parseProcessingTime(efile, verbose=0)
                print("generateLevelNext %s: legacy:  read processing time from %s: %.1f" % (tag, efile, dtx))
                researchOA.writeprocessingTime(pfile, dtx)
                # create numbers file
            print(
                "generateLevelNext %s: legacy: rename: %s -> %s"
                % (tag, os.path.join(outputdir, edir, mdfilelegacy), os.path.join(outputdir, edir, mdfileX0))
            )
            os.rename(os.path.join(outputdir, edir, mdfilelegacy), os.path.join(outputdir, edir, mdfileX0))
            researchOA.touchFile(legacyfile, "converted legacy format: %s\n" % tag)

        if verbose >= 1:
            print("generateLevelTwo: %s: already complete (old format)" % (tag))
        # raise Exception('old style top level directory: %s' % (join(outputdir, edir) ))

        joblist = []
        runcomplete = None
        cmd = None
        return joblist, runcomplete, cmd

    zipfile = join(outputdir, edir, "splitted-%s.zip" % splitTag(lvls))
    numbers = join(outputdir, edir, numbersFile(lvls))
    if checkFiles([numbers, zipfile]):
        if verbose >= 1:
            print("generateLevelNext: %s: already complete (zipfile)" % (tag))
        joblist = []
        runcomplete = None
        cmd = None
        return joblist, runcomplete, cmd

    joblist = []

    os.chdir(outputdir)
    if not researchOA.checkFilesOA(afilee):
        print("generateLevelNext: %s: note: file %s does not exist!" % (tag, afilee))
        runcomplete = None
        cmd = None
        return joblist, runcomplete, cmd
        pass  # raise Exception('file missing')

    myq = "q72h"

    if len(lvls) >= 2:
        myq = "q7d"

    """ Step 4a """
    if N > 56:
        qtime = 240
    else:
        qtime = 5
    print("# generateLevelNext: step 4a: tag %s, priorextend %s " % (tag, priorextend))
    if 1:
        edirf = edir

        lfile = os.path.join(edir, "lock-extend-prior.txt")
        # new style (sub-level 3)
        cmdlogfile = os.path.join(edir, rfile.replace(".oa", "-extend.txt"))
        maxk = splitdata["levels"][level + 1]
        checkstart = None
        cmd = ""
        if priorextend:
            cmd += "cd %s; mkdir -p %s;" % (outputdir, edirf)
            cmd += os.linesep + "touch %s" % lfile
            cmd += os.linesep + 'echo "discardJ5 %d" > %s\n' % (discardJ5, join(edirf, "extendoptions.txt"))
            if splitmode:
                cmd += (
                    os.linesep
                    + "oa_depth_extend --maxk %d -Q %d -v 2 -j 1 -m 9 --discardJ5 %d --split -f Z %s -o %s | tee %s;"
                    % (maxk, qtime, discardJ5, afilee, ebasefull, cmdlogfile)
                )
            else:
                cmd += (
                    os.linesep
                    + "oa_depth_extend --maxk %d -Q %d -v 2 -j 1 -m 9 --discardJ5 %d -f Z %s -o %s | tee %s;"
                    % (maxk, qtime, discardJ5, afilee, ebasefull, cmdlogfile)
                )
            cmd += os.linesep + "rm -f %s" % lfile + os.linesep
            checkstart = afilee
        cmd += "cd %s; mkdir -p %s;" % (outputdir, edirf)
        cmd += os.linesep + "# split data at level %d" % (level) + os.linesep
        splitcmd, splitfile, splitout, done = doSplitFile(
            lvls, splitdata, adata, verbose=verbose >= 2, outputdir=outputdir
        )
        cmd += splitcmd
        if checkstart is None:
            checkstart = splitfile
        splitcheckfile = join(splitdir(lvls), splitname(lvls + [splitdata[level]["n"] - 1]))
        if priorextend:
            shorttag = "ES%s"
            jobtype = "extend and split %s"
        else:
            shorttag = "S%s"
            jobtype = "split %s"
        j = job(
            cmd=cmd,
            jobtype=jobtype % tag,
            shorttag=shorttag % tag,
            ncores=1,
            checkfiles=[splitcheckfile],
            checkfilesstart=[checkstart],
        )
        joblist.append(j)

        if j.complete():
            if verbose:
                print("# generateLevelNext %s: extend and split complete" % (tag,))

            for kk in range(splitdata[level]["n"]):
                rfile = splitname(lvls + [kk])
                edir0X = splitdir(lvls)
                tagX = splitTag(lvls + [kk])
                edirX = splitdir(lvls + [kk])
                cmdlogfile = os.path.join(edirX, rfile.replace(".oa", "-extend.txt"))
                ebasefull = os.path.join(edirX, rfile).replace(".oa", "-extend")
                afilee = os.path.join(edir0X, rfile)
                # mdfileX='md5check-%s.txt' % tagX
                mdfileX = researchOA.mdFile(lvls + [kk])
                mdfileX2 = "md5check-%s.txt" % tagX.replace(".", "-")
                lfile = os.path.join(edirX, "lock-extend-%s.txt" % tagX)

                if os.path.exists(os.path.join(outputdir, lfile)):
                    print("lock file for %s exists, skipping" % tagX)
                    continue

                # cmd=''
                cmd = os.linesep + "# extend data at level (%d+?): %s" % (level, (tagX) + os.linesep)
                maxk = splitdata["levels"][level + 2]
                cmd += "cd %s; \nmkdir -p %s;\n" % (outputdir, edirX)
                cmd += os.linesep + 'echo "discardJ5 %d" > %s\n' % (discardJ5, join(edirX, "extendoptions.txt"))
                cmd += os.linesep + "touch %s" % lfile
                cmd += (
                    os.linesep
                    + "oa_depth_extend --maxk %d -Q %d -v 2 -j 1 -m 9 --discardJ5 %d --split %d -f Z %s -o %s | tee %s;"
                    % (maxk, qtime, discardJ5, splitmode, afilee, ebasefull, cmdlogfile)
                )
                cmd += os.linesep + "sleep .1"
                cmd += os.linesep + "rm -f %s" % lfile
                cmd += os.linesep + "cd %s; gzip -q -f *.oa; " % os.path.join(outputdir, edirX)
                cmd += os.linesep + "md5sum *.oa.gz > %s; " % (mdfileX,)
                j = job(
                    cmd=cmd,
                    checkfiles=[join(edirX, mdfileX)],
                    checkfilesstart=[afilee],
                    queue=myq,
                    ncores=ncores,
                    jobtype="E%s" % tagX,
                    shorttag="E%s" % tagX,
                )
                # fix
                if os.path.exists(join(edirX, mdfileX2)):
                    pass
                else:
                    if verbose >= 2:
                        print("# generateLevelNext %s: extend %d" % (tag, kk))
                    joblist.append(j)

                if 0:
                    print("xxx: adding job %s" % str(j))
                    print("  start %s" % (str(j.checkfilesstart)))
                    print("  check %s" % (str(j.checkfiles)))
                    print(os.getcwd())
                    exit(0)
        # for each splitted file do an extend
        # gather results
    runcomplete = None
    cmd = None
    return joblist, runcomplete, cmd


# %%


def makeJobList(scriptdir, jobs, ncores, queue, verbose=1):
    slist = []
    nj = 0
    for idx, j in enumerate(jobs):
        if not j.complete() and j.canrun():
            if ncores > 0:
                j.ncores = ncores
            outfile, substr = createJobScript(j, index=idx, verbose=0, queue=queue, scriptdir=scriptdir)
            slist += [substr]
            nj += 1
            if verbose >= 2:
                print("job: %s" % j.shorttag)
            if verbose >= 2 and idx < 5:
                print("  checkfiles: %s" % (j.checkfiles))
            if verbose >= 4:
                j.canrun(verbose=1)
                print("  checkfilesstart: %s" % (j.checkfilesstart))
        else:
            if verbose >= 3:
                print("job: %s: complete %d, canrun %d" % (j.shorttag, j.complete(), j.canrun()))
            if verbose >= 3 and idx < 1:
                j.canrun(verbose=1)
                if verbose >= 4:
                    print("  checkfilesstart: %s" % (j.checkfilesstart))
                print("  checkfiles: %s" % (j.checkfiles))

    print("created %d/%d jobs in %s" % (nj, len(jobs), scriptdir))

    fid = open(join(scriptdir, "subs.sh"), "wt")
    for i, s in enumerate(slist):
        _ = fid.write("%s\n" % s)
    fid.close()


# %%


def checkNumbers(basedir, lvls=[1, 95], verbose=1, debug=False):
    """Compare the results of a numbers file and statistics file

    Raises an error if the results are unequal
    """
    splitdir = researchOA.splitdir(lvls)
    checkdir = os.path.join(basedir, splitdir)  # ,/sp0-split-1/sp1-split-95'
    nfile = os.path.join(checkdir, researchOA.numbersFile(lvls))
    jfile = os.path.join(checkdir, researchOA.numbersFile(lvls, tag="jstats"))
    if debug:
        nfile = os.path.join(basedir, researchOA.numbersFile(lvls))
        jfile = os.path.join(basedir, researchOA.numbersFile(lvls, tag="jstats"))

    tag = splitTag(lvls)

    if verbose:
        print("checkNumbers: read %s" % nfile)

    try:
        mtime = os.path.getmtime(nfile)
    except OSError:
        mtime = 0
    xdate = datetime.datetime.fromtimestamp(mtime)
    xdate0 = datetime.datetime(2016, 1, 1)
    if xdate < xdate0 and mtime > 0:
        print("numbers file for %s is old: %s" % (lvls, time.ctime(mtime)))

    ndata = readNumbersFile(nfile)

    jdata = oapackage.readStatisticsFile(jfile, verbose=verbose >= 1)
    if verbose >= 3:
        print("--- jdata ")
        jdata.show()
    if not jdata.isOpen():
        raise Exception("file %s not found!" % jfile)

    tt = jdata.getTotals()

    for x in ndata:
        k = x[0]
        narrays = x[1]
        if verbose >= 2:
            print("checkNumbers %s: k %d: narrays %d" % (tag, k, narrays))

        if k >= len(tt) and narrays == 0 and k > 22:
            if verbose >= 2:
                print("high value of k, probably everying is zero, skipping")
            continue
        if tt[k] != narrays:
            raise Exception("number of arrays does not match! k %d: numbers %d, jstats %d" % (k, narrays, tt[k]))

    return ndata, jdata
