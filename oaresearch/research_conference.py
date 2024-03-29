""" Research topics related to conference designs

Pieter Eendebak <pieter.eendebak@gmail.com>

"""

import collections
import copy
import json
import operator
import os
import pickle
import shutil
import sys
import time
from collections import OrderedDict
from functools import lru_cache, reduce
from typing import Any, List, Tuple, Union

import json_tricks
import numpy as np
from joblib import Parallel, delayed

try:
    pass
except BaseException:
    pass

from collections import namedtuple

import oapackage
import oapackage.conference
import oapackage.markup as markup
import oapackage.oahelper as oahelper
from oapackage.markup import oneliner as e
from oapackage.oahelper import create_pareto_element

import oaresearch
import oaresearch.filetools
from oaresearch.misc_utils import index_sorted_array
from oaresearch.research import citation


#%%
class hashable_array:
    r"""Hashable wrapper for ndarray objects.

    Instances of ndarray are not hashable, meaning they cannot be added to
    sets, nor used as keys in dictionaries. This is by design - ndarray
    objects are mutable, and therefore cannot reliably implement the
    __hash__() method.

    The hashable class allows a way around this limitation. It implements
    the required methods for hashable objects in terms of an encapsulated
    ndarray object. This can be either a copied instance (which is safer)
    or the original object (which requires the user to be careful enough
    not to modify it).
    """

    def __init__(self, wrapped, tight=False):
        r"""Creates a new hashable object encapsulating an ndarray.

        wrapped
            The wrapped ndarray.

        tight
            Optional. If True, a copy of the input ndaray is created.
            Defaults to False.
        """
        self.__tight = tight
        self.__wrapped = np.array(wrapped) if tight else wrapped
        self.__hash = hash(self.__wrapped.tobytes())

    def __eq__(self, other):
        return np.all(self.__wrapped == other.__wrapped)

    def __hash__(self):
        return self.__hash

    def unwrap(self):
        r"""Returns the encapsulated ndarray.

        If the wrapper is "tight", a copy of the encapsulated ndarray is
        returned. Otherwise, the encapsulated ndarray itself is returned.
        """
        if self.__tight:
            return np.array(self.__wrapped)

        return self.__wrapped


def make_hashable_array(design):
    return hashable_array(np.array(design))


# %%

DesignType = Union[oapackage.array_link, np.ndarray]


def conference_design_extensions(
    array: DesignType, verbose: int = 0, filter_symmetry: bool = True, conference_generator=None
) -> List:
    """Return list of all extensions of a conference design

    All extensions are generated, minus symmetry conditions.
    For details see `oapackage.generateSingleConferenceExtensions`

    Args:
        arrays : List of arrays to be extended
        verbose: Verbosity level
        filter_symmetry: If True then discard any designs that are now minimal according to the row permutations
            from the row symmetry group of the input design
    Returns:
        List of extensions

    """
    j1zero = 0
    if not isinstance(array, oapackage.array_link):
        array = oapackage.makearraylink(array)
    conference_type = oapackage.conference_t(array.n_rows, array.n_columns, j1zero)

    zero_index = -1
    filterj2 = 1
    filterj3 = 0
    filter_symmetry = True  # we can use symmetry reduction, since any the other filtering is not related to the symmetry of the design

    if conference_generator is None:
        extensions = oapackage.generateSingleConferenceExtensions(
            array, conference_type, zero_index, verbose >= 2, filter_symmetry, filterj2, filterj3, filter_symmetry
        )
    else:
        # TODO assert proper conference class and symmetry types
        extensions = conference_generator.generateCandidates(array)
    extensions = [oapackage.hstack(array, extension) for extension in extensions]

    return list(extensions)


def conference_design_has_extension(design: DesignType) -> bool:
    """Return true if specified conference design has an extension with an additional column

    Args:
        design: Conference design
    Returns:
        True of the design has extensions
    """
    design_np = np.array(design)
    ee = conference_design_extensions(design_np)
    return len(ee) > 0


#%%

DesignStack = namedtuple("DesignStack", ["designs", "nauty_designs", "nauty_sort_index", "nauty_designs_sorted"])


def reduce_minimal_form(design, design_stack):
    """Reduce design to minimal form using full stack of designs"""
    all_data, all_data_nauty, sort_indices, nauty_designs_sorted = design_stack
    nauty_form = oapackage.reduceConference(design)
    k = nauty_form.shape[1]

    sort_idx = sort_indices[k]
    if nauty_designs_sorted is None:
        nauty_sorted = [all_data_nauty[k][idx] for idx in sort_idx]
    else:
        nauty_sorted = nauty_designs_sorted[k]
    sorted_idx = index_sorted_array(nauty_sorted, nauty_form)
    idx = sort_idx[sorted_idx]

    return all_data[k][idx]


@lru_cache(maxsize=400000)
def cached_conference_design_has_maximal_extension(design, verbose=0, Nmax=None, conference_generator=None):
    design_stack = conference_design_has_maximal_extension.design_stack
    if design_stack is None:
        raise Exception("initialize the design stack first!")

    if isinstance(design, hashable_array):
        design = design.unwrap()

    design_np = np.array(design)
    N = design_np.shape[0]
    if Nmax is None:
        Nmax = N
    k = design_np.shape[1]

    if not N == design_stack[0][4][0].shape[0]:
        raise Exception("N {N} does not match design stack")

    if k == Nmax:
        return True

    ee = conference_design_extensions(design_np, conference_generator=conference_generator)

    if verbose:
        print(f"design {N} {k}: check {len(ee)} extensions")
    result = False
    for subidx, extension_design in enumerate(ee):
        if verbose:
            print(f"design {N} {k}: subidx {subidx}")
        extension_design_link = oapackage.makearraylink(extension_design)
        md = reduce_minimal_form(extension_design_link, design_stack)

        if cached_conference_design_has_maximal_extension(
            make_hashable_array(md), verbose=verbose, Nmax=Nmax, conference_generator=conference_generator
        ):
            result = True
            break
    if verbose:
        print(f"design {N} {k}: {result}")
    return result


@oahelper.static_var("design_stack", None)
def conference_design_has_maximal_extension(design, verbose=0, Nmax=None, conference_generator=None) -> bool:
    """Determine whether a design has an extension to a maximal design

    The static variable design_stack needs to be initialized with a dictionary containing all designs
    with the number of rows specifed by the design.

    Args:
        design: Conference design
        N: Number of rows
    Returns:
        True if the design can be extended to the full number of columns
    """

    if not isinstance(design, hashable_array):
        design = make_hashable_array(design)
    result = cached_conference_design_has_maximal_extension(
        design, verbose=verbose, Nmax=Nmax, conference_generator=conference_generator
    )
    return result


def _flatten(data):
    if len(data) == 0:
        return data
    return reduce(operator.concat, data)


def _conference_design_extensions_numpy(array, verbose=0):
    lst = conference_design_extensions(array, verbose)
    return [np.array(array) for array in lst]


def extend_conference_designs_full(extensions: List, number_parallel_jobs=4) -> List:
    """Extend a list of conference designs with all possible extensions minus symmetry"""
    if number_parallel_jobs > 1:
        extensions_numpy = _flatten(
            Parallel(n_jobs=number_parallel_jobs)(
                delayed(_conference_design_extensions_numpy)(np.array(array)) for array in extensions
            )
        )
        extensions = [oapackage.makearraylink(array) for array in extensions_numpy]
    else:
        extensions = _flatten([conference_design_extensions(array) for array in extensions])
    return extensions


def maximal_extension_size(array: oapackage.array_link, verbose: int = 1) -> Tuple[int, list]:
    """Calculate maximum number of columns in an extention of specified design

    Args:
        design:
    Returns:
        Maximum number of columns
    """

    maximum_number_of_columns = array.n_columns

    N = array.n_rows

    extensions0 = [array]
    extensions = extensions0

    t0 = time.time()
    for c in range(extensions[0].n_columns, N):

        extensions2 = extend_conference_designs_full(extensions)
        if verbose >= 2:
            print(f"maximal_extension_size: N {N} columns {c}->{c+1}: {len(extensions)}->{len(extensions2)}")
        extensionsr = oapackage.selectConferenceIsomorpismClasses(extensions2, verbose=0)
        if verbose:
            dt = time.time() - t0
            print(
                f"maximal_extension_size: N {N} columns {c}->{c+1}: {len(extensions)}->{len(extensions2)}->{len(extensionsr)}, {dt:.1f} [s]"
            )
        if len(extensionsr) == 0:
            break

        extensions = extensionsr

        maximum_number_of_columns = c + 1

    return maximum_number_of_columns, extensions


# %%


def select_even_odd_conference_designs(cfile: str) -> Tuple:
    """Select the even-odd conference designs from a file with designs"""
    na = oapackage.nArrayFile(cfile)

    eolist: List[Any] = []
    if na > 100000:
        af = oapackage.arrayfile_t(cfile)
        for ii in range(na):
            if ii % (200 * 1e3) == 0 or ii == na - 1:
                print("select_even_odd_conference_designs: %s: %d/%d" % (cfile, ii, af.narrays))
            al = af.readnext()
            if ii == 0:
                if al.min() > -1:
                    raise Exception("not a conference matrix?!")
            if not oapackage.isConferenceFoldover(al):
                eolist += [al]
        af.closefile()
    else:
        ll = oapackage.readarrayfile(cfile)
        na = len(ll)
        if len(ll) > 0:
            if ll[0].min() > -1:
                raise Exception("not a conference matrix?!")

        eolist = [al for al in ll if not oapackage.isConferenceFoldover(al)]
    return na, eolist


# %%


def cdesignTag(
    N,
    kk,
    page,
    outputdir,
    tdstyle="",
    tags=["cdesign", "cdesign-diagonal", "cdesign-diagonal-r"],
    tagtype=["full", "r", "r"],
    verbose=1,
    ncache=None,
    subpage=None,
    generated_result=None,
    conference_html_dir=None,
):
    """Create html tag for oa page

    Args:
        N (int): number of rows
        kk (int): number of columns
        page (object):
        outputdir (str):
        tdstyle (str):
        tags (list):
        tagtype (list):
        verbose (int):
        ncache (dict): store results

    """
    cfile, nn, mode = conferenceResultsFile(N, kk, outputdir, tags=tags, tagtype=tagtype, verbose=1)

    if ncache is not None:
        if "full" not in ncache:
            ncache["full"] = {}
        ncache["full"]["N%dk%d" % (N, kk)] = nn

    cfilex = oapackage.oahelper.checkOAfile(cfile)
    if cfilex is not None:
        cfilebase = os.path.basename(cfilex)
    else:
        cfilebase = None

    if generated_result is not None and not len(generated_result) == 0:
        if generated_result["pareto_results"]["full_results"]:
            nn = generated_result["pareto_results"]["narrays"]
            cfile = None
            mode = "full"

    if page is not None:
        if subpage:
            hreflink = os.path.join("conference", subpage)
            if verbose:
                print("cdesignTag: hreflink: %s" % subpage)
        else:
            hreflink = "conference/%s" % cfilebase

        txt, hyper_link = htmlTag(nn, kk, N, mode=mode, href=hreflink, ncache=ncache, verbose=verbose >= 2)
        if verbose >= 2:
            print(
                "cdesignTag: N %d, k %d, html txt %s"
                % (
                    N,
                    kk,
                    txt,
                )
            )
        if hyper_link and (cfilex is not None) and not hreflink.endswith("html"):
            # no html page, just copy OA file
            if verbose >= 2:
                print("cdesignTag: N %d, ncols %d: copy OA file" % (N, kk))
            shutil.copyfile(cfilex, os.path.join(conference_html_dir, cfilebase))
        page.td(txt, style=tdstyle)
    else:
        if verbose >= 2:
            print(cfile)
    return cfile


def generate_even_odd_conference_designs(outputdir):
    """Select even-odd conference designs from generated double conference designs"""
    Nrange = range(0, 82, 2)  # full range
    # Nrange=range(44, 45, 2)
    # Nrange=range(4, 48, 2)
    Nrange = range(74, 82, 2)
    tag = "dconferencej1j3"
    for Ni, N in enumerate(Nrange):
        kmax = int(np.ceil(N / 2) + 1)
        krange = range(2, kmax)
        for ki, kk in enumerate(krange):

            cfile = cdesignTag(
                N,
                kk,
                page=None,
                outputdir=outputdir,
                tags=[tag, tag + "-r"],
                tagtype=["full", "r"],
                conference_html_dir=None,
            )

            na, eolist = select_even_odd_conference_designs(cfile)

            cfileout = cfile.replace(tag, tag + "-eo")
            print("  out: %s: %d -> %d" % (cfileout, na, len(eolist)))
            if 1:
                oapackage.writearrayfile(cfileout, eolist, oapackage.ABINARY_DIFF, N, kk)
                xfile = cfileout + ".gz"
                if os.path.exists(xfile):
                    print("removing file %s" % (xfile))
                    os.remove(xfile)
                if 1:
                    if len(eolist) > 100:
                        cmd = "gzip -f %s" % cfileout
                        os.system(cmd)


# %%


def reduce_single_conference(arrays, verbose=0):
    """Reduce a list of double conference arrays to single conference arrays

    Arrays that are not foldover arrays are discarded.

    Args:
        arrays (list): list of dobule conference designs
    Returns:
        list: list containing the corresponding single conference designs
    """
    narrays = len(arrays)
    arrays = [array for array in arrays if oapackage.isConferenceFoldover(array)]
    if verbose:
        print("reduce_single_conference: reduce %d arrays to %d single conference designs" % (narrays, len(arrays)))

    def reduce_single(array):
        Nsingle = int(array.n_rows / 2)
        perm = oapackage.double_conference_foldover_permutation(array)
        return oapackage.array_link(np.array(array)[perm[0:Nsingle], :])

    arrays = [reduce_single(array) for array in arrays]
    return arrays


class SingleConferenceParetoCombiner:
    def __init__(self, outputdir, cache_dir, cache=False, verbose=1, pareto_method_options=None):
        """Class to generate statistics and Pareto optimality results for a conference design class from double conference designs"""
        self.outputdir = outputdir
        self.cache_dir = cache_dir
        self.cache = cache
        self.verbose = verbose
        if pareto_method_options is None:
            pareto_method_options = {}
        self._pareto_method_options = pareto_method_options

    def append_basepath(self, afile):
        return os.path.join(self.outputdir, afile)

    def pareto_file(self, filename):
        pfile = os.path.join(self.cache_dir, filename)
        oapackage.mkdirc(os.path.split(pfile)[0])
        return pfile

    def stats_file(self, filename):
        pfile = os.path.join(self.cache_dir, filename).replace(".oa", ".json")
        oapackage.mkdirc(os.path.split(pfile)[0])
        return pfile

    def combined_results_file(self, number_columns):
        pfile = os.path.join(self.cache_dir, "combined-single-conference-pareto-results-k%d.json" % number_columns)
        oapackage.mkdirc(os.path.split(pfile)[0])
        return pfile

    def pre_calculate(self, arrayfiles):

        for ii, afile in enumerate(arrayfiles):
            outputfile = self.pareto_file(afile)
            outputfile_stats = self.stats_file(afile)

            if os.path.exists(outputfile) and self.cache:
                continue

            oapackage.oahelper.tprint("ParetoCalculator: pre_calculate file %d/%d: %s" % (ii, len(arrayfiles), afile))
            arrays = oapackage.readarrayfile(self.append_basepath(afile))
            number_arrays = len(arrays)
            arrays = reduce_single_conference(arrays, verbose=1)

            presults, _ = oaresearch.research_conference.calculateConferencePareto(
                arrays, **self._pareto_method_options
            )

            pareto_designs = [oapackage.array_link(array) for array in presults["pareto_designs"]]
            print("generate %s: %d arrays" % (outputfile, len(pareto_designs)))
            oapackage.writearrayfile(outputfile, oapackage.arraylist_t(pareto_designs), oapackage.ABINARY)
            with open(outputfile_stats, "w") as fid:
                json.dump({"number_arrays": number_arrays, "number_conference_arrays": len(arrays)}, fid)

    @staticmethod
    def combine_statistics(stats, extra_stats):
        if stats is None:
            return copy.copy(extra_stats)
        combined_stats = copy.copy(stats)
        for field in ["number_arrays", "number_conference_arrays"]:
            combined_stats[field] = stats[field] + extra_stats[field]

        return combined_stats

    def write_combined_results(self, number_columns, results):
        results["pareto_designs"] = [np.array(array) for array in results["pareto_designs"]]
        with open(self.combined_results_file(number_columns), "w") as fid:
            json_tricks.dump(results, fid, indent=4)

    def load_combined_results(self, number_columns):
        with open(self.combined_results_file(number_columns)) as fid:
            results = json_tricks.load(fid)
        return results

    def calculate(self, arrayfiles):
        """Calculate statistics over generated designs

        Args:
            lst (list): list of files with designs
        """

        pareto_arrays = []
        combined_stats = None
        for afile in arrayfiles:
            oapackage.oahelper.tprint("ParetoCalculator: calculate %s" % afile)
            outputfile = self.pareto_file(afile)
            outputfile_stats = self.stats_file(afile)

            arrays = oapackage.readarrayfile(outputfile)
            pareto_arrays += list(arrays)

            stats = json.load(open(outputfile_stats))
            combined_stats = self.combine_statistics(combined_stats, stats)

        presults, _ = oaresearch.research_conference.calculateConferencePareto(
            pareto_arrays, **self._pareto_method_options
        )
        # remove invalid fields
        for tag in ["B4", "F4"]:
            if tag + "_max" in presults:
                presults.pop(tag + "_max")
        for tag in ["rankinteraction", "ranksecondorder"]:
            if tag + "_min" in presults:
                presults.pop(tag + "_min")
        presults["combined_statistics"] = combined_stats

        presults["_pareto_method_options"] = self._pareto_method_options

        return presults


# %%


def generate_or_load_conference_results(
    N,
    number_of_columns,
    outputdir,
    dc_outputdir,
    double_conference_cases=(),
    addExtensions=True,
    addMaximumExtensionColumns=False,
):
    """Calculate results for conference designs class

    In data is either calculated directly, or loaded from pre-generated data gathered from double conference designs.

    """
    pareto_results = OrderedDict({"N": N, "ncolumns": number_of_columns, "full_results": 0, "no_results": True})

    from_double_conference = N in double_conference_cases

    addExtensions = N <= 26
    pareto_method_options = {
        "verbose": 1,
        "addProjectionStatistics": None,
        "addExtensions": addExtensions,
        "addMaximumExtensionColumns": addMaximumExtensionColumns,
    }

    if from_double_conference:
        if number_of_columns > N:
            return pareto_results, None

        print("generate_or_load_conference_results: N %d: loading from doubleconference results" % (N,))
        dc_dir = os.path.join(dc_outputdir, "doubleconference-%d" % (2 * N))
        if not os.path.exists(dc_dir):
            return {}, None
        cache_dir = oapackage.mkdirc(os.path.join(dc_dir, "sc_pareto_cache"))
        pareto_calculator = SingleConferenceParetoCombiner(
            dc_dir, cache_dir=cache_dir, cache=True, pareto_method_options=None
        )
        pareto_results = pareto_calculator.load_combined_results(number_of_columns)
        pareto_results["narrays"] = pareto_results["combined_statistics"]["number_conference_arrays"]
        pareto_results["idstr"] = "cdesign-%d-%d" % (pareto_results["N"], pareto_results["ncolumns"])

        pareto_results["full"] = True
        pareto_results["full_results"] = True
        pareto_results["_from_double_conference"] = True
        cfile = None
    else:
        cfile, nn, mode = conferenceResultsFile(
            N,
            number_of_columns,
            outputdir,
            tags=["cdesign", "cdesign-diagonal", "cdesign-diagonal-r"],
            tagtype=["full", "r", "r"],
            verbose=1,
        )
        ll = oapackage.readarrayfile(cfile)
        narrays = len(ll)

        if mode == "full" or narrays < 1000:
            presults, pareto = calculateConferencePareto(ll, N=N, k=number_of_columns, **pareto_method_options)
            pareto_results = generateConferenceResults(presults, ll, ct=None, full=(mode == "full"))
            pareto_results["arrayfile"] = cfile
        else:
            cfile = None
    return pareto_results, cfile


# %%


def generateConference(
    N, kmax=None, verbose=1, diagc=False, nmax=None, selectmethod="random", tag="cdesign", outputdir=None
):
    """Generate sequece of conference designs

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
        tag += "-diagonal"
    if nmax is not None:
        tag += "-r"

    al = ctype.create_root()

    ll = oapackage.arraylist_t()
    ll.push_back(al)
    LL = [[]] * (kmax)
    LL[1] = ll
    print("generateConference: start: %s" % ctype)
    if outputdir is not None:
        _ = oapackage.writearrayfile(os.path.join(outputdir, "cdesign-%d-%d.oa" % (N, 2)), LL[1], oapackage.ATEXT, N, 2)

    for extcol in range(2, kmax):
        if verbose:
            print("generateConference: extcol %d: %d designs" % (extcol, len(LL[extcol - 1])))
            sys.stdout.flush()
        LL[extcol] = oapackage.extend_conference(LL[extcol - 1], ctype, verbose=verbose >= 2)

        LL[extcol] = oapackage.selectConferenceIsomorpismClasses(LL[extcol], verbose >= 1)

        LL[extcol] = oapackage.sortLMC0(LL[extcol])

        if nmax is not None:
            na = min(nmax, len(LL[extcol]))
            if na > 0:
                if selectmethod == "random":
                    idx = np.random.choice(len(LL[extcol]), na, replace=False)
                    LL[extcol] = [LL[extcol][i] for i in idx]
                elif selectmethod == "first":
                    LL[extcol] = [LL[extcol][i] for i in range(na)]
                else:
                    # mixed case
                    raise Exception("not implemented")
        afmode = oapackage.ATEXT
        if len(LL[extcol]) > 1000:
            afmode = oapackage.ABINARY
        if outputdir is not None:
            _ = oapackage.writearrayfile(
                os.path.join(outputdir, "%s-%d-%d.oa" % (tag, N, extcol + 1)), LL[extcol], afmode, N, extcol + 1
            )

    ll = [len(l) for l in LL]
    if verbose:
        print("generated sequence: %s" % ll)
    return LL


# %%
def conferenceJ4(al, jj=4):
    """Calculate J4 values for a conference matrix"""

    al = oapackage.makearraylink(al)
    return oapackage.Jcharacteristics_conference(al, number_of_columns=jj)


@oapackage.oahelper.deprecated
def conferenceSecondOrder(al, include_so=False):
    """Calculate second-order interaction matrix for a conference matrix"""
    x = np.array(al)
    k = al.n_columns

    if include_so:
        offset = 0
        m = int(k * (k + 1) / 2)
    else:
        offset = 1
        m = int(k * (k - 1) / 2)
    y = np.zeros((x.shape[0], m))
    idx = 0
    for ii in range(k):
        for jj in range(ii + offset, k):
            y[:, idx] = x[:, ii] * x[:, jj]
            idx = idx + 1
    return y


def conferenceStatistics(al, verbose=0):
    """Calculate statistics for a conference design

    Args:
        al (array): design to use
    Returns:
        list: f4, b4, rank X2, rank X2 with quadratics
    """

    f4 = al.FvaluesConference(4)
    N = al.n_rows
    j4 = conferenceJ4(al)
    b4 = np.sum(np.array(j4) ** 2) / N**2

    ncols = al.n_columns
    modelmatrix = oapackage.array2modelmatrix(al, "i")[:, (1 + ncols) :]
    rank = np.linalg.matrix_rank(modelmatrix)

    modelmatrix_quadratic = oapackage.array2modelmatrix(al, "q")[:, (1 + ncols) :]
    rankq = np.linalg.matrix_rank(modelmatrix_quadratic)

    if verbose:
        print(f"f4: {f4}")
        print(f"j4: {j4}")
        print(f"rank X2: {rank}")
        print(f"rank X2+quadratics: {rankq}")
    return [f4, b4, rank, rankq]


def test_confJ4():
    al = oapackage.exampleArray(18)
    J = conferenceJ4(al)
    assert np.sum(np.abs(np.array(J)) == 12) == 1
    assert np.sum(np.abs(np.array(J)) == 0) == 23


def createConferenceParetoElement(
    al,
    addFoldover=True,
    addProjectionStatistics=True,
    addMaximality=False,
    addMaximumExtensionColumns=False,
    pareto=None,
    rounding_decimals=3,
):
    """Create Pareto element from conference design

    Returns:
        Tuple with pareto_element, data
    """
    rr = conferenceStatistics(al, verbose=0)
    [f4, b4, rankinteraction, ranksecondorder] = rr[0:4]
    f4minus = [-float(x) for x in f4]
    values = [[float(ranksecondorder)], [float(rankinteraction)], list(f4minus), [-float(b4)]]
    data = OrderedDict(ranksecondorder=ranksecondorder)
    data["rankinteraction"] = rankinteraction
    data["F4"] = f4
    data["B4"] = b4

    if addProjectionStatistics:
        proj_data = np.zeros((2, 3))
        proj_data[0] = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=4)
        proj_data[1] = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=5)

        proj_data = np.around(proj_data, rounding_decimals)

        for tag_index, tag in enumerate(["PEC", "PIC", "PPC"]):
            for ni, kk in enumerate([4, 5]):
                values += [[proj_data[ni, tag_index]]]

                data[tag + "%d" % kk] = proj_data[ni, tag_index]
    else:
        for tag in ["PEC", "PIC", "PPC"]:
            for kk in [4, 5]:
                data[tag + "%d" % kk] = None

    if addFoldover:
        foldover = oapackage.isConferenceFoldover(al)
        values += [[int(foldover)], [int(not foldover)]]
        data["foldover"] = int(foldover)
        data["notfoldover"] = int(not foldover)

    if addProjectionStatistics:
        assert len(values) == len(data.keys())
    if addMaximality:
        data["has_extensions"] = conference_design_has_extensions(al)
    if addMaximumExtensionColumns:
        data["maximum_extension_size"] = maximal_extension_size(al)[0]

    if pareto is None:
        pareto = oapackage.ParetoMultiDoubleLong()
    pareto_element = create_pareto_element(values, pareto=pareto)

    return pareto_element, data


@oapackage.oahelper.deprecated
def makePareto(presults, addFoldover=True):
    pareto = oapackage.ParetoMultiDoubleLong()

    for ii in range(len(presults.ranks)):
        f4minus = tuple(-x for x in presults.f4s[ii])
        values = [[int(presults.ranks[ii])], list(f4minus), [-presults.b4s[ii]]]
        if addFoldover:
            values += [[int(presults.foldover[ii])], [int(not presults.foldover[ii])]]
        val = create_pareto_element(values, pareto=pareto)

        pareto.addvalue(val, ii)

    return pareto


class pareto_results_structure(collections.OrderedDict):
    """Class to hold results of Pareto calculations"""

    def add_value(self, tag, value):
        mintag = tag + "_min"
        maxtag = tag + "_max"
        if mintag not in self:
            self[mintag] = value
        if maxtag not in self:
            self[maxtag] = value
        self[mintag] = min(value, self[mintag])
        self[maxtag] = max(value, self[maxtag])


def conference_design_has_extensions(array, verbose=0):
    """Return True if a single conference design has extensions"""
    j1zero = 0
    conference_type = oapackage.conference_t(array.n_rows, array.n_columns, j1zero)

    zero_index = -1
    filterj2 = 1
    filterj3 = 0
    filter_symmetry = (
        1  # we can use symmetry reduction, since any the other filtering is not related to the symmetry of the design
    )
    extensions = oapackage.generateSingleConferenceExtensions(
        array, conference_type, zero_index, verbose >= 2, filter_symmetry, filterj2, filterj3, filter_symmetry
    )

    result = len(extensions) > 0
    if verbose:
        print("conference_design_has_extensions: %s, found %d extensions" % (result, len(extensions)))
        if verbose >= 2:
            oapackage.showCandidates(extensions)

    return result


def conferenceParetoIdentifier():
    return "0.4"


def calculateConferencePareto(
    ll,
    N=None,
    k=None,
    verbose=1,
    add_data=True,
    addProjectionStatistics=None,
    addExtensions=False,
    addMaximumExtensionColumns=False,
    number_parallel_jobs=1,
):
    """Calculate Pareto optimal designs from a list of designs

    Args:
        ll (list): list of designs
        N (int)
        k (int)
        verbose (int)
        add_data (bool)
        addProjectionStatistics (None or bool)
        addExtensions (bool)
    Returns:
        presults, pareto
    """
    t0 = time.time()

    if verbose:
        print(
            "calculateConferencePareto: analysing %d arrays, addProjectionStatistics %s, addExtensions %s"
            % (len(ll), addProjectionStatistics, addExtensions)
        )

    if len(ll) > 0:
        N = ll[0].n_rows
        k = ll[0].n_columns

    presults = pareto_results_structure({"pareto_designs": []})
    pareto = oapackage.ParetoMultiDoubleLong()
    if N is None:
        presults["N"] = None
        return presults, pareto

    if addProjectionStatistics is None:
        if N <= 20:
            addProjectionStatistics = True
        else:
            addProjectionStatistics = False

    data = None
    t0 = time.time()

    def pareto_worker(al):
        al = oapackage.makearraylink(al)
        pareto_element, data = createConferenceParetoElement(
            al, addFoldover=False, addProjectionStatistics=addProjectionStatistics, pareto=None
        )
        pareto_element = [list(e.values) for e in list(pareto_element)]
        return pareto_element, data

    block_size = 200
    blocks = [(ii, min(len(ll), ii + block_size)) for ii in np.arange(0, len(ll), block_size)]

    def add_extra_data(presults, data, addProjectionStatistics):
        if add_data:
            for tag in ["ranksecondorder", "rankinteraction", "B4", "F4"]:
                presults.add_value(tag, data[tag])
            if addProjectionStatistics:
                for tag in ["PEC4", "PIC4"]:
                    presults.add_value(tag, data[tag])

    if number_parallel_jobs > 1:
        for blockidx, block in enumerate(blocks):
            dt = time.time() - t0
            oapackage.oahelper.tprint(
                f"calculateConferencePareto: N {N} column {k}: block {blockidx}/{len(blocks)}: ({dt:.1f} [s]): {str(pareto).strip()}"
            )

            xx = Parallel(n_jobs=number_parallel_jobs)(
                [delayed(pareto_worker)(np.array(al)) for al in ll[block[0] : block[1]]]
            )

            for jj, al in enumerate(ll[block[0] : block[1]]):
                pareto_element, data = xx[jj]

                pareto_element = oapackage.vector_mvalue_t_double(
                    [oapackage.mvalue_t_double(x) for x in pareto_element]
                )

                ii = int(jj + block[0])
                pareto.addvalue(pareto_element, ii)
                add_extra_data(presults, data, addProjectionStatistics)

    else:
        for ii, al in enumerate(ll):
            oapackage.oahelper.tprint(
                "calculateConferencePareto: N %s column %s: array %d/%d (%.1f [s]): %s"
                % (str(N), str(k), ii, len(ll), time.time() - t0, str(pareto).strip()),
                dt=2,
            )
            pareto_element, data = createConferenceParetoElement(
                al, addFoldover=False, addProjectionStatistics=addProjectionStatistics, pareto=None
            )

            pareto.addvalue(pareto_element, ii)
            add_extra_data(presults, data, addProjectionStatistics)

    presults["N"] = N
    presults["ncolumns"] = k

    if len(ll) > 0:
        presults.N = ll[0].n_rows
        presults.ncolumns = ll[0].n_columns

    if data is None:
        presults["pareto_type"] = "no design"
    else:
        presults["pareto_type"] = ", ".join([key for key in data.keys() if data[key] is not None])

    presults["pareto_type"] = presults["pareto_type"].replace("ranksecondorder", "r(2FI, QE)")
    presults["pareto_type"] = presults["pareto_type"].replace("rankinteraction", "r(2FI)")

    pareto.show()

    presults["pareto_indices"] = pareto.allindices()
    presults["nclasses"] = pareto.number()
    presults["npareto"] = pareto.numberindices()
    presults["_version"] = conferenceParetoIdentifier()

    presults["pareto_designs"] = [ll[ii] for ii in presults["pareto_indices"]]
    presults["pareto_data"] = []

    def pareto_worker(al):
        al = oapackage.makearraylink(al)

        pareto_element, data = createConferenceParetoElement(
            al, addFoldover=True, addMaximality=addExtensions, addMaximumExtensionColumns=addMaximumExtensionColumns
        )
        return data

    if number_parallel_jobs > 1:
        data_results = Parallel(n_jobs=number_parallel_jobs)(
            [delayed(pareto_worker)(np.array(al)) for al in presults["pareto_designs"]]
        )

        for ii, data in enumerate(data_results):
            presults["pareto_data"].append(data)

    else:
        for ii, al in enumerate(presults["pareto_designs"]):
            data = pareto_worker(al)
            presults["pareto_data"].append(data)

    presults["_pareto_processing_time"] = time.time() - t0
    presults = OrderedDict(presults)
    return presults, pareto


def test_calculateConferencePareto():
    ll = [oapackage.exampleArray(idx) for idx in [45, 46, 47, 45]]
    presults, _ = calculateConferencePareto(ll, N=None, k=None, verbose=1, add_data=True)


if __name__ == "__main__":
    test_calculateConferencePareto()


def showMaxZ(LL):
    """For a list of generated designs show the maximum zero position"""
    N = LL[3][0].n_rows

    for ii, L in enumerate(LL):
        k = ii + 1
        s = [oapackage.maxz(al) for al in L]
        mm, _ = np.histogram(s, range(N + 1))
        print("%d cols: maxz seq %s" % (k, list(mm)))


def generate_conference_latex_tables(htmlsubdir, verbose=1):
    """Generate LaTeX results tables from pre-generated result files"""
    for N in range(8, 25, 2):
        lst = oapackage.findfiles(htmlsubdir, "conference-N%d.*pickle" % N)
        if verbose:
            print("latex table: N %d: %d files" % (N, len(lst)))
        table = None

        kk = [oapackage.scanf.sscanf(file, "conference-N%dk%d")[1] for file in lst]
        lst = [lst[idx] for idx in np.argsort(kk)]

        for file in lst:
            r = pickle.load(open(os.path.join(htmlsubdir, file), "rb"))

            ncolumns = r["ncolumns"]
            rtable = r["rtable"]
            if rtable.size == 0:
                continue
            column = np.vstack((["k"], ncolumns * np.ones((rtable.shape[0] - 1, 1), dtype=int)))
            rtable = np.hstack((column, rtable))
            if table is None:
                table = rtable
            else:
                rtable = rtable[1:]
                table = np.vstack((table, rtable))
            # r['ncolumns']
        print(table)
        if len(lst) == 0:
            print("no results for N=%d" % N)
            continue

        offset_columns = [1, 2]
        for row in range(1, table.shape[0]):
            for col in offset_columns:
                table[row, col] = str(int(table[row, col]) + 1)
        latextable = oapackage.array2latex(
            table, hlines=[0], comment=["conference desgins N=%d" % (N), "offset for indices is 1"]
        )
        if verbose:
            print(latextable)
        with open(os.path.join(htmlsubdir, "conference-N%d-overview.tex" % (N,)), "w") as fid:
            fid.write(latextable)


def conferenceResultsFile(
    N, kk, outputdir, tags=["cdesign", "cdesign-diagonal", "cdesign-diagonal-r"], tagtype=["full", "r", "r"], verbose=1
):
    """Create html tag for oa page

    Args:
        N (int): number of rows
        kk (int): number of columns
        outputdir (str):
        tags (list):
        tagtype (list):
        verbose (int):
        ncache (dict): store results
    """
    for ii, tag in enumerate(tags):

        cfile0 = "%s-%d-%d.oa" % (tag, N, kk)
        cfile = os.path.join(outputdir, cfile0)
        gfile = os.path.join(outputdir, cfile0 + ".gz")
        if verbose >= 2:
            print("cdesignTag: try file %s" % cfile0)
        if os.path.exists(os.path.join(outputdir, cfile0)) and os.path.exists(gfile):
            oapackage.nArrayFile(cfile)
            oapackage.nArrayFile(gfile)
            raise Exception("both .oa and .oa.gz exist: %s" % cfile)

        nn = oapackage.nArrays(cfile)
        mode = tagtype[ii]
        cfilex = oapackage.oahelper.checkOAfile(cfile)
        if cfilex is not None:
            os.path.basename(cfilex)
        else:
            pass
        if nn >= 0:
            break

    if verbose:
        print("cdesignTag: N %d, kk %d: selected tag %s: nn %d" % (N, kk, tag, nn))
    # special case
    if kk == N and tag == "cdesign-diagonal":
        mode = "full"

    if verbose >= 2:
        print(cfile)
    return cfile, nn, mode


def generateConferenceResults(presults, ll, ct=None, full=None):
    pareto_results = presults
    pareto_results["type"] = "conference designs"
    pareto_results["arrayfile"] = None
    pareto_results["presults"] = None

    pareto_results["full"] = full
    pareto_results["full_results"] = full
    pareto_results["idstr"] = "cdesign-%d-%d" % (pareto_results["N"], pareto_results["ncolumns"])
    if ct is not None:
        pareto_results["ctidstr"] = ct.idstr()
        assert ct.N == pareto_results["N"]
        assert ct.ncols == pareto_results["ncolumns"]

    pareto_results["narrays"] = len(ll)

    pareto_results["pareto_designs"] = [np.array(array) for array in presults["pareto_designs"]]
    return pareto_results


# %% Webpage generation


def nprevzero(N, k, ncache):
    """Return true if any previous result was zero"""
    for ix in range(k - 1, 2, -1):
        # print(ix)
        p = ncache["full"].get("N%dk%d" % (N, ix), -1)
        # p = ncache['full'].get('N%dk%d' % (N, ix), -1)
        if p == 0:
            return True
    return False


def htmlTag(nn, kk, N, mode="full", href=None, ncache=None, verbose=0):
    """Create html tag for number of designs

    Args:
        nn (int): number of arrays
        kk (int): number of columns
        N (int): number of rows
        mode (str)
        href (None or str): hyperlink to subpage

    Returns:
        txt (str): link text
        hyper_link (bool): True if the txt is a hyperlink
    """
    hyper_link = False
    if nn >= 0:
        if mode == "full":
            txt = "%d" % nn
        else:
            if nn == 0:
                txt = "?"
            else:
                txt = "&ge; %d" % nn
        if href is not None:
            if nn < 6000 and nn > 0 or (href.endswith("html")):
                ss = e.a(txt, href=href, style="text-decoration: none;")
                hyper_link = True
                txt = ss
        else:
            pass
    else:
        if verbose:
            print("htmlTag: nn is negative")
        if kk <= N:
            if ncache is None:
                txt = "?"
            else:
                if verbose >= 1:
                    print("htmlTag: mode %s, N %d, k %d" % (mode, N, kk))
                if nprevzero(N, kk, ncache):
                    if verbose >= 1:
                        print("htmlTag: nprevzero(%d, %d, ncache) is True" % (N, kk))
                    txt = ""
                else:
                    txt = "?"
        else:
            txt = ""
    return txt, hyper_link


def latexResults(outputdir):
    X = []
    print("make latex results table...")

    NN = range(4, 31, 2)
    kk = range(0, np.max(NN) + 1)

    X = np.zeros((1 + 1 + len(kk) - 2, 1 + len(NN)), dtype=object)
    X[:] = ""
    X[0, 0] = ""
    X[1, 0] = "$k$"
    for ii, N in enumerate(NN):
        X[1, 1 + ii] = N
        X[0, 1 + ii] = ""
        for ki, k in enumerate(range(2, N + 1)):
            if k > N:
                X[1 + 1 + ki, 1 + ii] = ""
            else:
                cfile0 = "cdesign-%d-%d.oa" % (N, k)
                nn = oapackage.nArrays(os.path.join(outputdir, cfile0))
                if nn < 0 and k == N:
                    cfile0 = "cdesign-diagonal-%d-%d.oa" % (N, k)
                    nn = oapackage.nArrays(os.path.join(outputdir, cfile0))
                if nn < 0:
                    cfile0 = "cdesign-diagonal-%d-%d.oa" % (N, k)
                    nnm = oapackage.nArrays(os.path.join(outputdir, cfile0))
                    if nnm > 0:
                        X[1 + 1 + ki, 1 + ii] = r"$\ge %d$" % nnm
                    else:
                        X[1 + 1 + ki, 1 + ii] = "?"
                else:
                    X[1 + 1 + ki, 1 + ii] = nn
                X[1 + 1 + ki, 0] = "%d" % k
    X[0, 1] = "$N$"

    X[2, 0] = X[2, 0] + r"\rule{0pt}{2.9ex}"
    return X


def createConferenceDesignsPageHeader(page, makeheader, conference_class, ncolumns, full_results=False):
    xstr = "C(%d, %d)" % (conference_class.N, ncolumns)
    xstrplain = xstr
    if makeheader:
        page.init(
            title="Class %s" % xstrplain,
            css=("../oastyle.css"),
            lang="en",
            htmlattrs=dict({"xmlns": "http://www.w3.org/1999/xhtml", "xml:lang": "en"}),
            header="<!-- Start of page -->",
            htmlheader=oaresearch.research.oaCssStyle(addframe=True),
            bodyattrs=dict({"style": "padding-left: 3px;"}),
            doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
            metainfo=(
                {
                    "text/html": "charset=utf-8",
                    "keywords": "conference designs",
                    "robots": "index, follow",
                    "description": "conference designs",
                }
            ),
            footer="<!-- End of page -->",
        )

    if full_results:
        page.h1("Conference designs %s " % xstr)
    else:
        page.h1("Conference designs %s (<b>partial results</b>) " % xstr)
    oap = e.a("Orthogonal Array package", href=r"http://www.pietereendebak.nl/oapackage/index.html")
    pstr = "This page contains information about conference designs. "
    pstr += "The results have been generated with the %s." % oap
    pstr += " If you use these results, please cite the paper " + citation("conference", style="full") + "."
    page.p(pstr)


def createConferenceDesignsPageResultsTable(page, pareto_results, verbose=0):
    full_results = pareto_results.get("full")
    if full_results:
        page.h2("Results")
    else:
        page.h2("Results (partial results)")

    page.table()
    page.tr(style="font-weight: bold;")
    page.td("Statistic", style="padding-right:30px;")
    page.td(("Results"), style="padding-right:8px;")
    page.tr.close()

    def simpleRow(a, b):
        page.tr(style="")
        page.td(a, style="padding-right:30px;")
        page.td(b, style="padding-right:8px;")
        page.tr.close()

    narrays = pareto_results["narrays"]
    simpleRow("Number of non-isomorphic designs", str(pareto_results["narrays"]))

    if narrays > 0:
        if "ranksecondorder_min" in pareto_results:
            simpleRow(
                "Minimum/Maximum rank of model matrix with 2FI and QE",
                "%d/%d" % (pareto_results["ranksecondorder_min"], pareto_results["ranksecondorder_max"]),
            )
            simpleRow(
                "Minimum/Maximum rank of model matrix for 2FI",
                "%d/%d" % (pareto_results["rankinteraction_min"], pareto_results["rankinteraction_max"]),
            )
        else:
            simpleRow("Maximum rank of model matrix with 2FI and QE", "%d" % (pareto_results["ranksecondorder_max"]))
            simpleRow("Maximum rank of model matrix for 2FI", "%d" % (pareto_results["rankinteraction_max"]))
        if "B4_min" in pareto_results:
            simpleRow("Minimum B4", "%.4f" % pareto_results["B4_min"])
        if "B4_max" in pareto_results:
            simpleRow("Maximum B4", "%.4f" % pareto_results["B4_max"])
        if "F4_min" in pareto_results:
            simpleRow("Minimum F4", "{}".format(pareto_results["F4_min"]))
        if "F4_max" in pareto_results:
            simpleRow("Maximum F4", "{}".format(pareto_results["F4_max"]))
    if "totaltime" in list(pareto_results.keys()):
        simpleRow("Processing time", str(pareto_results["totaltime"]) + "s")

    if pareto_results.get("datafile_tag", None) is not None:
        simpleRow("Data", pareto_results["datafile_tag"])

    page.table.close()
    page.p(style="font-size: smaller;")
    page.add("Note: 2FI: two-factor interaction; QE: quadratic effect")
    page.p.close()


def createConferenceDesignsPageLoadDesignsFile(pareto_results, htmlsubdir=None, verbose=1):
    havearrayfile = 0
    if "arrayfile" in list(pareto_results.keys()):
        if pareto_results["arrayfile"] is not None:
            havearrayfile = 1

    if verbose:
        print("createConferenceDesignsPageLoadDesignsFile: havearrayfile %d" % havearrayfile)
    if havearrayfile:
        iafile = pareto_results["arrayfile"]
        outfile0 = pareto_results["idstr"] + ".oa"

        na = oapackage.nArrayFile(iafile)
        if verbose >= 2:
            print("conferenceDesignsPage: read arrayfile %s: na %d" % (iafile, na))

        if na < 5000 and na >= 0:
            if htmlsubdir is None:
                raise Exception("need html subdirectory to copy .oa file")
            outfilefinal = oaresearch.filetools.copyOAfile(
                iafile, htmlsubdir, outfile0, convert="T", zipfile=None, verbose=1, cache=0
            )

            if pareto_results.get("full", False):
                if verbose:
                    print("conferenceDesignsPage: full results")
                pareto_results["datafilestr"] = "all arrays"
                htag = oaresearch.htmltools.formatArrayHyperlink(pareto_results["datafilestr"], outfile0, iafile)

                pareto_results["datafilestr"] = "all arrays"
                pareto_results["datafile_tag"] = htag
            else:
                na = oapackage.nArrayFile(os.path.join(htmlsubdir, outfilefinal))
                pareto_results["datafilestr"] = "%d array(s)" % na
                htag = oaresearch.htmltools.formatArrayHyperlink(pareto_results["datafilestr"], outfile0, iafile)
                pareto_results["datafile_tag"] = htag

        else:
            if verbose:
                print("conferenceDesignsPage: no datafile (na %d)" % na)
            pareto_results["datafilestr"] = "-"
            pareto_results["datafile_tag"] = None


def createConferenceDesignsPageParetoTable(page, pareto_results, verbose=0, htmlsubdir=None):
    """Create table with Pareto results and add to the markup object

    Args:
        page (markup.page): html page to add table to

    Returns:
        rtable (array): generated table
    """
    if verbose:
        print("createConferenceDesignsPageParetoTable: start")

    pareto_indices = pareto_results["pareto_indices"]
    pareto_data = pareto_results["pareto_data"]

    add_extension_information = False
    add_maximum_extension_size = False

    if len(pareto_data) > 0:
        if pareto_data[0].get("has_extensions", None) is not None:
            add_extension_information = True
        if pareto_data[0].get("maximum_extension_size", None) is not None:
            add_maximum_extension_size = True

    if pareto_results["narrays"] > 0 and pareto_results.get("full_results"):
        add_extra = True
        if verbose:
            print("createConferenceDesignsPageParetoTable: do statistics2htmltable")

        header = ["Index Pareto file", "Index design file", "Rank 2FI and QE", "Rank 2FI", "F4", "B4"]
        if add_extra:
            for tag in ["PEC", "PIC", "PPC"]:
                for kk in [4, 5]:
                    header += [tag + "%d" % kk]
        if add_extension_information:
            header += ["Extensions"]
        if add_maximum_extension_size:
            header += ["Max. columns"]

        rtable = np.zeros((1 + len(pareto_results["pareto_indices"]), len(header)), dtype="|U208")
        rtable[:] = " "
        for ii, h in enumerate(header):
            rtable[0, ii] = header[ii]

        sort_indices = oapackage.sortrows(np.array([p["F4"] for p in pareto_results["pareto_data"]]))

        for ii, sort_index in enumerate(sort_indices):
            pareto_idx = sort_index
            array_list_idx = pareto_indices[sort_index]
            rank_secondorder = str(pareto_data[pareto_idx]["ranksecondorder"])
            rank_interaction = str(pareto_data[pareto_idx]["rankinteraction"])
            rowdata = [
                "%d" % pareto_idx,
                "%d" % array_list_idx,
                rank_secondorder,
                rank_interaction,
                str(pareto_data[pareto_idx]["F4"]),
                "%.2f" % (pareto_data[pareto_idx]["B4"]),
            ]
            rtable[ii + 1, 0 : len(rowdata)] = rowdata
            column_offset = len(rowdata)
            if add_extra:
                for tag in ["PEC", "PIC", "PPC"]:
                    for kk in [4, 5]:
                        rtable[ii + 1, column_offset] = "%.3f" % (pareto_data[pareto_idx][tag + "%d" % kk])
                        column_offset = column_offset + 1
            if add_extension_information:
                rtable[ii + 1, column_offset] = "Yes" if pareto_data[pareto_idx]["has_extensions"] > 0 else "No"
                column_offset = column_offset + 1
            if add_maximum_extension_size:
                rtable[ii + 1, column_offset] = pareto_data[pareto_idx]["maximum_extension_size"]
                column_offset = column_offset + 1

        subpage = oaresearch.research.array2html(
            rtable,
            header=1,
            tablestyle="border-collapse: collapse;",
            trclass="",
            tdstyle="padding-right:1em;",
            trstyle="",
            thstyle="text-align:left; padding-right: 1em;",
            comment=None,
        )
        page.br(clear="both")
        page.h2("Pareto optimal designs")
        page.p()
        if pareto_results["nclasses"] == 1:
            pareto_classes_text = "in %d class" % pareto_results["nclasses"]
        else:
            pareto_classes_text = "in %d classes" % pareto_results["nclasses"]
        if pareto_results["npareto"] == 1:
            page.add("There is %d Pareto optimal design %s." % (pareto_results["npareto"], pareto_classes_text))
        else:
            page.add("There are %d Pareto optimal designs %s." % (pareto_results["npareto"], pareto_classes_text))
        pareto_type = pareto_results["pareto_type"]
        if "," in pareto_type:
            k = pareto_type.rfind(", ")
            pareto_type = pareto_type[:k] + ", and " + pareto_type[k + 1 :]

        page.add("Pareto optimality is according to %s (any other statistics are ignored)." % pareto_type)
        page.p.close()
        if pareto_results.get("pareto_designs", None) is not None:
            pdesigns = pareto_results.get("pareto_designs", None)

        pfile0 = pareto_results["idstr"] + "-pareto.oa"

        if htmlsubdir is not None:
            pfile = os.path.join(htmlsubdir, pfile0)
            oapackage.writearrayfile(pfile, [oapackage.array_link(array) for array in pdesigns])
            page.p("All %s" % e.a("Pareto optimal designs", href=pfile0) + ".")

        page.add(str(subpage))
        page.br(clear="both")
    else:
        rtable = np.zeros((0, 0))

    return rtable


def _convert_to_latex_table(rtable, N, ncolumns, offset_columns=[0, 1]):
    import copy

    rtable_latex = copy.deepcopy(rtable)
    if len(rtable) > 1:
        for row in range(1, rtable_latex.shape[0]):
            for col in offset_columns:
                rtable_latex[row, col] = str(int(rtable_latex[row, col]) + 1)
    latextable = oapackage.array2latex(
        rtable_latex, hlines=[0], comment=["conference desgins N=%d, k=%d" % (N, ncolumns), "offset for indices is 1"]
    )
    return latextable


def conferenceDesignsPage(
    pareto_results, verbose=1, makeheader=True, htmlsubdir=None, generate_latex=True, html_template=False
):
    """Generate html page for class conference designs

    Args:
        pareto_results (dict): structure with results for Pareto optimal designs
        html_template (bool): If True then place html files in subfolder templates/
    """

    N = pareto_results["N"]
    ncolumns = pareto_results["ncolumns"]

    if makeheader:
        pass
    else:
        pass

    conference_class = oapackage.conference_t(N, ncolumns, 0)

    if verbose:
        print("conferenceDesignsPage: start of generation")

    page = markup.page()
    createConferenceDesignsPageHeader(
        page,
        makeheader=makeheader,
        conference_class=conference_class,
        ncolumns=ncolumns,
        full_results=pareto_results["full_results"],
    )
    createConferenceDesignsPageLoadDesignsFile(pareto_results, htmlsubdir=htmlsubdir)
    createConferenceDesignsPageResultsTable(page, pareto_results, verbose=verbose)
    rtable = createConferenceDesignsPageParetoTable(page, pareto_results, verbose=verbose, htmlsubdir=htmlsubdir)

    latextable = _convert_to_latex_table(rtable, N, ncolumns)

    if generate_latex:
        with open(os.path.join(htmlsubdir, "conference-N%dk%d.tex" % (N, ncolumns)), "w") as fid:
            fid.write(latextable)
        import pickle

        with open(os.path.join(htmlsubdir, "conference-N%dk%d-rtable.pickle" % (N, ncolumns)), "wb") as fid:
            pickle.dump({"rtable": rtable, "N": pareto_results["N"], "ncolumns": pareto_results["ncolumns"]}, fid)

    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)

    page.p("<br/>\n")
    page.p("Page generated on %s." % dstr)

    pstr = str(page).replace(
        '<meta content="charset=utf-8" name="text/html" />',
        '<meta http-equiv="Content-Type" content="charset=utf-8" name="text/html" />',
    )

    return pstr
