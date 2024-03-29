"""

Example script for calculating extension properties of conference designs

Pieter Eendebak <pieter.eendebak@gmail.com>
"""


# %% Load necessary packages
import os
import pickle
import platform
from os.path import join
from typing import Dict, Tuple

import numpy as np
import oapackage
import oapackage.graphtools

import oaresearch.research_conference
from oaresearch.research_conference import (
    conferenceDesignsPage,
    generate_or_load_conference_results,
)

generate_webpage = True


# %% Setup directories
resultsdir = join(os.path.expanduser("~"), "oatmp")

# outputdir of the code generating the designs
outputdir = oapackage.mkdirc(os.path.join(os.path.expanduser("~"), "oatmp", "conf"))

if platform.node() == "woelmuis":
    paperdir = "/home/eendebakpt/misc/oa/article-conference/"

if generate_webpage:
    htmldir = os.path.join(os.path.expanduser("~"), "misc", "oapage2")
    html_template = True

    oapackage.mkdirc(os.path.join(htmldir))
    conference_html_dir = oapackage.mkdirc(os.path.join(htmldir, "conference"))


def conference_designs_result_file(outputdir: str, tag: str, N: int, kk: int) -> str:
    cache_tag = f"result-pareto-v1-{oaresearch.research_conference.conferenceParetoIdentifier()}"
    oapackage.mkdirc(os.path.join(outputdir, cache_tag))
    cachefile = os.path.join(outputdir, cache_tag, tag + "-" + "N%dk%d" % (N, kk) + ".pickle")
    return cachefile


analyse_full = False
verbose = 1

# %% Generate subpages for the designs


def conferenceSubPages(
    tag="conference",
    Nmax=40,
    Nstart=4,
    kmax=None,
    outputdir=None,
    conference_html_dir=None,
    verbose=1,
    specials={},
    Nstep=2,
    NmaxPareto=40,
    cache=True,
    double_conference_cases=(24,),
    html_template=False,
    addMaximumExtensionColumns=False,
):
    """Generate subpages for single conference results

    Args:
        page (markup object)
        tag (str): type of data
        Nmax, Nstart, kmax (int)
        ta (str): text-align property
        verbose (int): verbosity level

    """

    debugdata = {}

    Nrange = range(Nstart, Nmax + 1, Nstep)
    if kmax == -2:
        kmax = int(np.ceil(Nmax / 2) + 1)
    if kmax is None:
        kmax = np.max(Nrange) + 2
    krange = range(2, kmax)
    if verbose:
        print("conferenceSubPages: tag %s, Nmax %d, kmax %d" % (tag, Nmax, kmax))

    subpages = {}
    subpages[tag] = {}

    for ki, kk in enumerate(krange):
        for Ni, N in enumerate(Nrange):

            if tag == "conference":
                if kk > N:
                    continue

            subpages[tag]["N%dk%d" % (N, kk)] = {}
            cachefile = conference_designs_result_file(outputdir, tag, N, kk)

            if cache and os.path.exists(cachefile):
                with open(cachefile, "rb") as fid:
                    if verbose >= 2:
                        print("loading results from cachefile %s" % cachefile)
                    pareto_results, cfile = pickle.load(fid)
                if verbose:
                    print("conferenceSubPages %s: from cache: N %d, columns %d" % (tag, N, kk))
                if pareto_results["presults"] is not None:
                    pversion = pareto_results["_version"]
                    cversion = oaresearch.research_conference.conferenceParetoIdentifier()
                    if pareto_results["_version"] != cversion:
                        raise Exception(
                            f"conference Pareto definition was updated! file {cachefile} data {pversion} code {cversion}"
                        )

            else:
                pareto_results, cfile = generate_or_load_conference_results(
                    N,
                    kk,
                    outputdir,
                    dc_outputdir=resultsdir,
                    double_conference_cases=double_conference_cases,
                    addMaximumExtensionColumns=N <= 20,
                )

                if verbose >= 2:
                    print("storing results in cachefile %s" % cachefile)

                with open(cachefile, "wb") as fid:
                    pareto_results["presults"] = None
                    pickle.dump((pareto_results, cfile), fid)

            if pareto_results.get("narrays", None) is None:
                continue

            # create HTML page
            page = conferenceDesignsPage(
                pareto_results,
                verbose=verbose >= 2,
                makeheader=True,
                htmlsubdir=conference_html_dir,
                html_template=html_template,
            )

            # write results
            htmlfile0 = tag + "-%d-%d.html" % (N, kk)
            if html_template:
                oapackage.mkdirc(os.path.join(conference_html_dir, "templates"))
                htmlfile = os.path.join(conference_html_dir, "templates", htmlfile0 + ".template")
            else:
                htmlfile = os.path.join(conference_html_dir, htmlfile0)

            sx = subpages[tag]["N%dk%d" % (N, kk)]
            sx["htmlpage0"] = htmlfile0
            sx["htmlpage"] = htmlfile
            sx["pareto_results"] = pareto_results
            sx["arrayfile"] = cfile

            if verbose >= 2:
                print("writing to %s" % htmlfile)
            with open(htmlfile, "w") as fid:
                fid.write(page)

    conferenceSubPages._debugdata = debugdata
    return subpages


generated_subpages = conferenceSubPages(
    tag="cdesign",
    Nmax=22,
    Nstart=4,
    verbose=1,
    cache=True,
    outputdir=outputdir,
    conference_html_dir=conference_html_dir,
    html_template=html_template,
)

#%%
from oaresearch.research_conference import (
    DesignStack,
    cached_conference_design_has_maximal_extension,
    conference_design_extensions,
    conference_design_has_maximal_extension,
    make_hashable_array,
)


def load_designs(N, k):
    data = generated_subpages["cdesign"][f"N{N}k{k}"]

    pr = data["pareto_results"]
    data["arrayfile"]

    cfile = os.path.join(outputdir, pr["idstr"] + ".oa")
    designs = oapackage.readarrayfile(cfile)
    return designs


def load_design_stack(Nx: int) -> Tuple[Dict]:
    """Load all conference designs for a specified number of rows"""
    all_data = {}
    all_data_nauty = {}
    sort_idx = {}
    data_nauty_sorted = {}
    for k in range(2, Nx + 1):
        designs = load_designs(Nx, k)

        all_data[k] = designs
        all_data_nauty[k] = [oapackage.reduceConference(d) for d in all_data[k]]

        sort_idx[k] = np.argsort(np.array(all_data_nauty[k], dtype=object))
        data_nauty_sorted[k] = [all_data_nauty[k][idx] for idx in sort_idx[k]]
    return DesignStack(all_data, all_data_nauty, sort_idx, data_nauty_sorted)


N = 16
all_data, all_data_nauty, sort_idx, _ = design_stack = load_design_stack(N)
k = 8

designs = load_designs(N, k)


from oaresearch.research_conference import conference_design_has_extension as has_extension

#%%
design = designs[10]  # .showarray()
ee = conference_design_extensions(design)

for d in ee:
    d.showarray()

#%%
design = designs[10]  # .showarray()
design = design.selectLastColumns(7)
ee = conference_design_extensions(design)

design.showarray()
print("extensions")
for d in ee:
    d.showarray()


#%%
conference_design_has_maximal_extension.design_stack = design_stack
cached_conference_design_has_maximal_extension.cache_clear()
print("first...")
conference_design_has_maximal_extension((designs[1]), verbose=1)
print("second...")
conference_design_has_maximal_extension((designs[1]), verbose=1)

#%%
t0 = oapackage.get_time_ms()
rr = []
design_has_extension = []
for idx, design in enumerate(designs):
    r = conference_design_has_maximal_extension(design, verbose=0)
    rr.append(r)
    design_has_extension.append(has_extension(design))
    print(f"{idx}: {r}")
dt = oapackage.get_time_ms() - t0

na = len(rr)
print(
    f"total time: {dt:.1f} [s]: N {N} k {k} max extension {np.sum(rr)}/{na}, has extension {np.sum(design_has_extension)}/{na}"
)

print(cached_conference_design_has_maximal_extension.cache_info())

#%%
for N in range(4, 20, 2):
    t0 = oapackage.get_time_ms()
    ctype = oapackage.conference_t(N, N, 0)
    conference_generator = oapackage.CandidateGeneratorConference(oapackage.array_link(), ctype)
    # conference_generator=None
    design_stack = load_design_stack(N)
    conference_design_has_maximal_extension.design_stack = design_stack
    cached_conference_design_has_maximal_extension.cache_clear()

    for k in range(2, N):
        designs = load_designs(N, k)
        design_has_max_extension_results = []
        design_has_extension_results = []
        for idx, design in enumerate(designs):
            if verbose >= 2:
                print(f"N {N} k {k} idx {idx}")
            hd = make_hashable_array(design)
            design_has_max_extension = conference_design_has_maximal_extension(
                hd, verbose=0, conference_generator=conference_generator
            )
            design_has_max_extension_results.append(design_has_max_extension)
            if 0:
                if design_has_max_extension:
                    design_has_extension_results.append(True)
                else:
                    design_has_extension_results.append(has_extension(design))
        dt = oapackage.get_time_ms() - t0

        na = len(design_has_max_extension_results)
        print(
            f"total time: {dt:.1f} [s]: N {N} k {k}: # designs with maximal extension {np.sum(design_has_max_extension_results)}/{na}, # designs with an extension {np.sum(design_has_extension_results)}/{na}"
        )

#%%
N = 20
ctype = oapackage.conference_t(N, N, 0)
conference_generator = oapackage.CandidateGeneratorConference(oapackage.array_link(), ctype)


design_stack = load_design_stack(N)
conference_design_has_maximal_extension.design_stack = design_stack

t0 = oapackage.get_time_ms()
k_values = [2, 3, 4, 5] + list(range(12, 20))
if analyse_full:
    k_values = list(range(2, 20))
for k in k_values:
    designs = load_designs(N, k)
    design_has_max_extension_results = []
    design_has_extension_results = []
    for idx, design in enumerate(designs):
        # print(f'N {N} k {k} idx {idx}')
        hd = make_hashable_array(design)
        design_has_max_extension = conference_design_has_maximal_extension(
            hd, verbose=0, conference_generator=conference_generator
        )
        design_has_max_extension_results.append(design_has_max_extension)
        if 0:
            if design_has_max_extension:
                design_has_extension_results.append(True)
            else:
                design_has_extension_results.append(has_extension(design))
    dt = oapackage.get_time_ms() - t0

    na = len(design_has_max_extension_results)
    print(
        f"total time: {dt:.1f} [s]: N {N} k {k}: # designs with maximal extension {np.sum(design_has_max_extension_results)}/{na}"
    )  # , # designs with an extension {np.sum(design_has_extension_results)}/{na}')

#%%
if analyse_full:
    N = 22
    design_stack = load_design_stack(N)
    conference_design_has_maximal_extension.design_stack = design_stack
    cached_conference_design_has_maximal_extension.cache_clear()

    ctype = oapackage.conference_t(N, N, 0)
    conference_generator = oapackage.CandidateGeneratorConference(oapackage.array_link(), ctype)

    kk = [10, 11]
    kk = range(2, 12)
    kk = kk[::-1]
    t0 = oapackage.get_time_ms()
    for k in kk:
        designs = load_designs(N, k)
        design_has_max_extension_results = []
        design_has_extension_results = []
        for idx, design in enumerate(designs):
            if verbose >= 2:
                if idx % 5000 == 0:
                    dt = oapackage.get_time_ms() - t0
                    print(f"N {N} k {k} idx {idx}/{len(designs)} dt {dt:.1f} [s]")
            hd = make_hashable_array(design)
            design_has_max_extension = conference_design_has_maximal_extension(
                hd, verbose=0, Nmax=12, conference_generator=conference_generator
            )
            design_has_max_extension_results.append(design_has_max_extension)
        dt = oapackage.get_time_ms() - t0

        na = len(design_has_max_extension_results)
        print(
            f"total time: {dt:.1f} [s]: N {N} k {k}: # designs with maximal extension {np.sum(design_has_max_extension_results)}/{na}"
        )


#%% N=22 pareto designs
N = 22
design_stack = load_design_stack(N)
conference_design_has_maximal_extension.design_stack = design_stack
cached_conference_design_has_maximal_extension.cache_clear()

ctype = oapackage.conference_t(N, N, 0)
conference_generator = oapackage.CandidateGeneratorConference(oapackage.array_link(), ctype)

#%%
basedir = "/home/eendebakpt/Dropbox/conference designs/designs/confpage_dc/conference/"
for k in range(7, 12 + 1):
    cfile = os.path.join(basedir, f"cdesign-{N}-{k}-pareto.oa")
    # cfile=os.path.join(basedir, f'cdesign-{N}-{k}.oa')
    pp = oapackage.readarrayfile(cfile)

    for idx, design in enumerate(pp):
        hd = make_hashable_array(design)
        design_has_max_extension = conference_design_has_maximal_extension(
            hd, verbose=0, Nmax=12, conference_generator=conference_generator
        )
        print(f"N {N} k {k} pareto design {idx}: design_has_max_extension {design_has_max_extension}")

print(cached_conference_design_has_maximal_extension.cache_info())

#%%

for idx, design in enumerate(designs):
    ee = conference_design_extensions(design)
    xx = [len(conference_design_extensions(d)) for d in ee]
    print(xx)

    r = len(ee)
    print(f"{idx}: {r}")

#%%

design_stack = load_design_stack(16)
designs = design_stack[0][6]
design = designs[0]
conference_design_has_maximal_extension.design_stack = design_stack
r = conference_design_has_maximal_extension(make_hashable_array(design), verbose=1)
print(r)

ee = conference_design_extensions(design)
xx = [len(conference_design_extensions(d)) for d in ee]
print(xx)


# ee[0].showarraycompact()

#%%
#%timeit has_maximal_extension(make_hashable_array(design), verbose=1)

#%% Determine number of extensible designs
from oaresearch.research_conference import conference_design_extensions

for Nx in range(4, 16, 2):
    for k in range(2, Nx + 1):
        data = generated_subpages["cdesign"][f"N{Nx}k{k}"]

        data.keys()

        pr = data["pareto_results"]
        data["arrayfile"]

        cfile = os.path.join(outputdir, pr["idstr"] + ".oa")
        # print(cfile)

        designs = oapackage.readarrayfile(cfile)

        ee = [len(conference_design_extensions(design)) > 0 for design in designs]
        # print(ee)
        ne = np.sum(ee)
        print(f"N {Nx} k {k}: {len(designs)} designs, {ne} can be extended")
