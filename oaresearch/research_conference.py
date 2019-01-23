# -*- coding: utf-8 -*-
""" Module to generate D-optimal designs

For more information see: https://doi.org/10.1080/00401706.2016.1142903

Pieter Eendebak <pieter.eendebak@gmail.com>

"""

from __future__ import print_function

import collections
import itertools
import os
import sys
import time
import pickle

import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    pass

import oapackage
import oapackage.conference
import oapackage.markup as markup
import oapackage.oahelper as oahelper
from oapackage.markup import oneliner as e

import oaresearch
from oaresearch.research import citation
import oaresearch.filetools

from oapackage.conference import momentMatrix, modelStatistics, conferenceProjectionStatistics

# %%


import json
import copy
import json_tricks


def select_even_odd_conference_designs(cfile):
    """ Select the even-odd conference designs from a file with designs """
    na = oapackage.nArrayFile(cfile)

    eolist = []
    if na > 100000:
        af = oapackage.arrayfile_t(cfile)
        for ii in range(na):
            if ii % (200 * 1e3) == 0 or ii == na - 1:
                print('select_even_odd_conference_designs: %d/%d' % (cfile, ii, af.narrays))
            al = af.readnext()
            if ii == 0:
                if al.min() > -1:
                    raise Exception('not a conference matrix?!')
            if not oapackage.isConferenceFoldover(al):
                eolist += [al]
        af.closefile()
    else:
        ll = oapackage.readarrayfile(cfile)
        na = len(ll)
        if len(ll) > 0:
            if ll[0].min() > -1:
                raise Exception('not a conference matrix?!')

        eolist = [al for al in ll if not oapackage.isConferenceFoldover(al)]
    return na, eolist


def reduce_single_conference(arrays, verbose=0):
    """ Reduce a list of double conference arrays to single conference arrays

    Arrays that are not foldover arrays are discarded.

    Args:
        arrays (list): list of dobule conference designs
    Returns:
        list: list containing the corresponding single conference designs
    """
    narrays = len(arrays)
    arrays = [array for array in arrays if oapackage.isConferenceFoldover(array)]
    if verbose:
        print('reduce_single_conference: reduce %d arrays to %d single conference designs' % (narrays, len(arrays)))

    def reduce_single(array):
        Nsingle = int(array.n_rows / 2)
        perm = oapackage.double_conference_foldover_permutation(array)
        return oapackage.array_link(np.array(array)[perm[0:Nsingle], :])

    arrays = [reduce_single(array) for array in arrays]
    return arrays


class SingleConferenceParetoCombiner:

    def __init__(self, outputdir, cache_dir, cache=False, verbose=1):
        """ Class to generate statistics and Pareto optimality results for a conference design class from double conference designs """
        self.outputdir = outputdir
        self.cache_dir = cache_dir
        self.cache = cache
        self.verbose = verbose

    def append_basepath(self, afile):
        return os.path.join(self.outputdir, afile)

    def pareto_file(self, filename):
        pfile = os.path.join(self.cache_dir, filename)
        oapackage.mkdirc(os.path.split(pfile)[0])
        return pfile

    def stats_file(self, filename):
        pfile = os.path.join(self.cache_dir, filename).replace('.oa', '.json')
        oapackage.mkdirc(os.path.split(pfile)[0])
        return pfile

    def combined_results_file(self, number_columns):
        pfile = os.path.join(self.cache_dir, 'combined-single-conference-pareto-results-k%d.json' % number_columns)
        oapackage.mkdirc(os.path.split(pfile)[0])
        return pfile

    def pre_calculate(self, arrayfiles):

        for ii, afile in enumerate(arrayfiles):
            outputfile = self.pareto_file(afile)
            outputfile_stats = self.stats_file(afile)

            if os.path.exists(outputfile) and self.cache:
                continue

            oapackage.oahelper.tprint('ParetoCalculator: pre_calculate file %d/%d: %s' % (ii, len(arrayfiles), afile))
            arrays = oapackage.readarrayfile(self.append_basepath(afile))
            number_arrays = len(arrays)
            arrays = reduce_single_conference(arrays, verbose=1)

            presults, pareto = oaresearch.research_conference.calculateConferencePareto(arrays)

            pareto_designs = [oapackage.array_link(array) for array in presults['pareto_designs']]
            print('generate %s: %d arrays' % (outputfile, len(pareto_designs)))
            oapackage.writearrayfile(outputfile, oapackage.arraylist_t(pareto_designs), oapackage.ABINARY)
            with open(outputfile_stats, 'wt') as fid:
                json.dump({'number_arrays': number_arrays, 'number_conference_arrays': len(arrays)}, fid)

    def combine_statistics(self, stats, extra_stats):
        if stats is None:
            return copy.copy(extra_stats)
        combined_stats = copy.copy(stats)
        for field in ['number_arrays', 'number_conference_arrays']:
            combined_stats[field] = stats[field] + extra_stats[field]

        return combined_stats

    def write_combined_results(self, number_columns, results):
        results['pareto_designs'] = [np.array(array) for array in results['pareto_designs']]
        with open(self.combined_results_file(number_columns), 'wt') as fid:
            json_tricks.dump(results, fid, indent=4)

    def load_combined_results(self, number_columns):
        with open(self.combined_results_file(number_columns), 'rt') as fid:
            results = json_tricks.load(fid)
        return results

    def calculate(self, arrayfiles):
        """ Calculate statistics over generated designs

        Args:
            lst (list): list of files with designs
        """

        pareto_arrays = []
        combined_stats = None
        for afile in arrayfiles:
            oapackage.oahelper.tprint('ParetoCalculator: calculate %s' % afile)
            outputfile = self.pareto_file(afile)
            outputfile_stats = self.stats_file(afile)

            arrays = oapackage.readarrayfile(outputfile)
            pareto_arrays += list(arrays)

            stats = json.load(open(outputfile_stats, 'rt'))
            combined_stats = self.combine_statistics(combined_stats, stats)

        presults, pareto = oaresearch.research_conference.calculateConferencePareto(pareto_arrays)
        # remove invalid fields
        for tag in ['B4', 'F4']:
            if tag + '_max' in presults:
                presults.pop(tag + '_max')
        for tag in ['rankinteraction', 'ranksecondorder']:
            if tag + '_min' in presults:
                presults.pop(tag + '_min')
        presults['combined_statistics'] = combined_stats
        return presults


# %%


def generate_or_load_conference_results(N, number_of_columns, outputdir, dc_outputdir, double_conference_cases=[]):
    """ Calculate results for conference designs class

    In data is either calculated directly, or loaded from pre-generated data gathered from double conference designs.

    """
    pareto_results = OrderedDict({'N': N, 'ncolumns': number_of_columns, 'full_results': 0, 'no_results': True})

    from_double_conference = N in double_conference_cases
    if from_double_conference:
        if number_of_columns > N:
            return pareto_results, None

        print('generate_or_load_conference_results: N %d: loading from doubleconference results' % (N,))
        dc_dir = os.path.join(dc_outputdir, 'doubleconference-%d' % (2 * N))
        cache_dir = oapackage.mkdirc(os.path.join(dc_dir, 'sc_pareto_cache'))
        pareto_calculator = SingleConferenceParetoCombiner(dc_dir, cache_dir=cache_dir, cache=True)
        pareto_results = pareto_calculator.load_combined_results(number_of_columns)
        pareto_results['narrays'] = pareto_results['combined_statistics']['number_conference_arrays']
        pareto_results['idstr'] = 'cdesign-%d-%d' % (pareto_results['N'], pareto_results['ncolumns'])

        pareto_results['full'] = True
        pareto_results['full_results'] = True
        pareto_results['_from_double_conference'] = True
        cfile = None
    else:
        cfile, nn, mode = conferenceResultsFile(N, number_of_columns, outputdir,
                                                tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
                                                tagtype=['full', 'r', 'r'], verbose=1)
        ll = oapackage.readarrayfile(cfile)
        narrays = len(ll)

        if mode == 'full' or narrays < 1000:
            presults, pareto = calculateConferencePareto(ll, N=N, k=number_of_columns, verbose=1,
                                                         addProjectionStatistics=None)
            pareto_results = generateConferenceResults(presults, ll, ct=None, full=(mode == 'full'))
            pareto_results['arrayfile'] = cfile
        else:
            cfile = None
    return pareto_results, cfile


# %%


def generateConference(N, kmax=None, verbose=1, diagc=False, nmax=None, selectmethod='random', tag='cdesign',
                       outputdir=None):
    """ Generate sequece of conference designs

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
    if nmax is not None:
        tag += '-r'

    al = ctype.create_root()

    ll = oapackage.arraylist_t()
    ll.push_back(al)
    LL = [[]] * (kmax)
    LL[1] = ll
    print('generateConference: start: %s' % ctype)
    if outputdir is not None:
        _ = oapackage.writearrayfile(
            os.path.join(outputdir, 'cdesign-%d-%d.oa' % (N, 2)), LL[1], oapackage.ATEXT, N, 2)

    for extcol in range(2, kmax):
        if verbose:
            print('generateConference: extcol %d: %d designs' % (extcol, len(LL[extcol - 1])))
            sys.stdout.flush()
        LL[extcol] = oapackage.extend_conference(
            LL[extcol - 1], ctype, verbose=verbose >= 2)

        LL[extcol] = oapackage.selectConferenceIsomorpismClasses(LL[extcol], verbose >= 1)

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
                    # mixed case
                    raise Exception('not implemented')
        afmode = oapackage.ATEXT
        if (len(LL[extcol]) > 1000):
            afmode = oapackage.ABINARY
        if outputdir is not None:
            _ = oapackage.writearrayfile(os.path.join(
                outputdir, '%s-%d-%d.oa' % (tag, N, extcol + 1)), LL[extcol], afmode, N, extcol + 1)

    ll = [len(l) for l in LL]
    if verbose:
        print('generated sequence: %s' % ll)
    return LL


# %%
def conferenceJ4(al, jj=4):
    """ Calculate J4 values for a conference matrix """

    al = oapackage.makearraylink(al)
    return oapackage.Jcharacteristics_conference(al, jj=jj)

    al = np.array(al)
    nc = al.shape[1]
    jj = []
    for c in list(itertools.combinations(range(nc), 4)):
        X = al[:, c]
        j4 = int(np.sum(np.prod(X, axis=1)))
        jj += [j4]
    return jj


@oapackage.oahelper.deprecated
def conferenceSecondOrder(al, include_so=False):
    """ Calculate second-order interaction matrix for a conference matrix """
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
    """ Calculate statistics for a conference design

    Args:
        al (array): design to use
    Returns:
        list: f4, b4, rank X2, rank X2 with quadratics
    """

    f4 = al.FvaluesConference(4)
    N = al.n_rows
    j4 = conferenceJ4(al)
    b4 = np.sum(np.array(j4) ** 2) / N ** 2

    ncols = al.n_columns
    modelmatrix = oapackage.array2modelmatrix(al, 'i')[:, (1 + ncols):]
    rank = np.linalg.matrix_rank(modelmatrix)

    modelmatrix_quadratic = oapackage.array2modelmatrix(al, 'q')[:, (1 + ncols):]
    rankq = np.linalg.matrix_rank(modelmatrix_quadratic)

    if verbose:
        print('f4: %s' % (f4,))
        print('j4: %s' % (j4,))
        print('rank X2: %s' % (rank,))
        print('rank X2+quadratics: %s' % (rankq,))
    return [f4, b4, rank, rankq]


def test_confJ4():
    al = oapackage.exampleArray(18)
    J = conferenceJ4(al)
    assert (np.sum(np.abs(np.array(J)) == 12) == 1)
    assert (np.sum(np.abs(np.array(J)) == 0) == 23)


from oapackage.oahelper import create_pareto_element
from collections import OrderedDict


def createConferenceParetoElement(al, addFoldover=True, addProjectionStatistics=True, pareto=None, rounding_decimals=3):
    """ Create Pareto element from conference design """
    rr = conferenceStatistics(al, verbose=0)
    [f4, b4, rankinteraction, ranksecondorder] = rr[0:4]
    f4minus = [-float(x) for x in f4]
    values = [[float(ranksecondorder)], [float(rankinteraction)], list(f4minus), [-float(b4)]]
    data = OrderedDict(ranksecondorder=ranksecondorder)
    data['rankinteraction'] = rankinteraction
    data['F4'] = f4
    data['B4'] = b4

    if addProjectionStatistics:
        proj_data = np.zeros((2, 3))
        proj_data[0] = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=4)
        proj_data[1] = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=5)

        proj_data = np.around(proj_data, rounding_decimals)

        for tag_index, tag in enumerate(['PEC', 'PIC', 'PPC']):
            for ni, kk in enumerate([4, 5]):
                values += [[proj_data[ni, tag_index]]]

                data[tag + '%d' % kk] = proj_data[ni, tag_index]
    else:
        for tag in ['PEC', 'PIC', 'PPC']:
            for kk in [4, 5]:
                data[tag + '%d' % kk] = None

    if addFoldover:
        foldover = oapackage.isConferenceFoldover(al)
        values += [[int(foldover)], [int(not foldover)]]
        data['foldover'] = int(foldover)
        data['notfoldover'] = int(not foldover)

    if addProjectionStatistics:
        assert (len(values) == len(data.keys()))
    pareto_element = create_pareto_element(values, pareto=pareto)

    return pareto_element, data


@oapackage.oahelper.deprecated
def makePareto(presults, addFoldover=True):
    pareto = oapackage.ParetoMultiDoubleLong()

    for ii in range(len(presults.ranks)):
        f4minus = tuple([-x for x in presults.f4s[ii]])
        values = [[int(presults.ranks[ii])], list(f4minus), [-presults.b4s[ii]]]
        if addFoldover:
            values += [[int(presults.foldover[ii])], [int(not presults.foldover[ii])]]
        val = create_pareto_element(values, pareto=pareto)

        pareto.addvalue(val, ii)

    return pareto


class pareto_results_structure(collections.OrderedDict):
    """ Class to hold results of Pareto calculations """

    def add_value(self, tag, value):
        mintag = tag + '_min'
        maxtag = tag + '_max'
        if not mintag in self:
            self[mintag] = value
        if not maxtag in self:
            self[maxtag] = value
        self[mintag] = min(value, self[mintag])
        self[maxtag] = max(value, self[maxtag])


def calculateConferencePareto(ll, N=None, k=None, verbose=1, add_data=True, addProjectionStatistics=None):
    """ Calculate Pareto optimal designs from a list of designs """
    if verbose:
        print('calculateConferencePareto: analysing %d arrays, addProjectionStatistics %s' % (
            len(ll), addProjectionStatistics))

    if len(ll) > 0:
        N = ll[0].n_rows
        k = ll[0].n_columns

    presults = pareto_results_structure({'pareto_designs': []})
    pareto = oapackage.ParetoMultiDoubleLong()
    if N is None:
        presults['N'] = None
        # pareto['nclasses=']=0
        return presults, pareto

    if addProjectionStatistics is None:
        if (N <= 20):
            addProjectionStatistics = True
        else:
            addProjectionStatistics = False

    data = None
    t0 = time.time()
    for ii, al in enumerate(ll):
        oapackage.oahelper.tprint('calculateConferencePareto: N %s column %s: array %d/%d (%.1f [s]): %s' %
                                  (str(N), str(k), ii, len(ll), time.time() - t0, str(pareto).strip()), dt=2)
        pareto_element, data = createConferenceParetoElement(al, addFoldover=False,
                                                             addProjectionStatistics=addProjectionStatistics,
                                                             pareto=pareto)

        pareto.addvalue(pareto_element, ii)

        if add_data:
            for tag in ['ranksecondorder', 'rankinteraction', 'B4', 'F4']:
                presults.add_value(tag, data[tag])
            if addProjectionStatistics:
                for tag in ['PEC4', 'PIC4']:
                    presults.add_value(tag, data[tag])

    presults['N'] = N
    presults['ncolumns'] = k

    if len(ll) > 0:
        presults.N = ll[0].n_rows
        presults.ncolumns = ll[0].n_columns

    if data is None:
        presults['pareto_type'] = 'no design'
    else:
        presults['pareto_type'] = ', '.join([key for key in data.keys() if data[key] is not None])
    if 0:
        if addProjectionStatistics:
            for kk in [4, 5]:
                presults['pareto_type'] += ', '
                presults['pareto_type'] += ', '.join(['PIC%d' % kk, 'PEC%d' % kk, 'PPC%d' % kk])

    presults['pareto_type'] = presults['pareto_type'].replace('ranksecondorder', 'r(2FI, QE)')
    presults['pareto_type'] = presults['pareto_type'].replace('rankinteraction', 'r(2FI)')

    pareto.show()

    presults['pareto_indices'] = pareto.allindices()
    presults['nclasses'] = pareto.number()
    presults['npareto'] = pareto.numberindices()
    presults['_version'] = '0.2'

    presults['pareto_designs'] = [ll[ii] for ii in presults['pareto_indices']]
    presults['pareto_data'] = []
    for ii, al in enumerate(presults['pareto_designs']):
        pareto_element, data = createConferenceParetoElement(al, addFoldover=True)
        presults['pareto_data'].append(data)

    presults = OrderedDict(presults)
    return presults, pareto


def test_calculateConferencePareto():
    ll = [oapackage.exampleArray(idx) for idx in [45, 46, 47, 45]]
    presults, _ = calculateConferencePareto(ll, N=None, k=None, verbose=1, add_data=True)


if __name__ == '__main__':
    test_calculateConferencePareto()


def showMaxZ(LL):
    """ For a list of generated designs show the maximum zero position """
    N = LL[3][0].n_rows

    for ii, L in enumerate(LL):
        k = ii + 1
        s = [oapackage.maxz(al) for al in L]
        mm, _ = np.histogram(s, range(N + 1))
        print('%d cols: maxz seq %s' % (k, list(mm)))


def generate_conference_latex_tables(htmlsubdir, verbose=1):
    """ Generate LaTeX results tables from pre-generated result files """
    for N in range(8, 25, 2):
        lst = oapackage.findfiles(htmlsubdir, 'conference-N%d.*pickle' % N)
        if verbose:
            print('latex table: N %d: %d files' % (N, len(lst)))
        table = None

        kk = [oapackage.scanf.sscanf(file, 'conference-N%dk%d')[1] for file in lst]
        lst = [lst[idx] for idx in np.argsort(kk)]

        for file in (lst):
            r = pickle.load(open(os.path.join(htmlsubdir, file), 'rb'))

            ncolumns = r['ncolumns']
            rtable = r['rtable']
            if rtable.size == 0:
                continue
            column = np.vstack((['k'], ncolumns * np.ones((rtable.shape[0] - 1, 1), dtype=int)))
            rtable = np.hstack((column, rtable))
            if table is None:
                table = rtable
            else:
                rtable = rtable[1:]
                table = np.vstack((table, rtable))
            # r['ncolumns']
        print(table)
        if len(lst) == 0:
            print('no results for N=%d' % N)
            continue

        offset_columns = [1, 2]
        for row in range(1, table.shape[0]):
            for col in offset_columns:
                table[row, col] = str(int(table[row, col]) + 1)
        latextable = oapackage.array2latex(table, hlines=[0],
                                           comment=['conference desgins N=%d' % (N), 'offset for indices is 1'])
        if verbose:
            print(latextable)
        with open(os.path.join(htmlsubdir, 'conference-N%d-overview.tex' % (N,)), 'wt') as fid:
            fid.write(latextable)


def conferenceResultsFile(N, kk, outputdir, tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
                          tagtype=['full', 'r', 'r'], verbose=1):
    """ Create html tag for oa page

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

        cfile0 = '%s-%d-%d.oa' % (tag, N, kk)
        cfile = os.path.join(outputdir, cfile0)
        gfile = os.path.join(outputdir, cfile0 + '.gz')
        if verbose >= 2:
            print('cdesignTag: try file %s' % cfile0)
        if os.path.exists(os.path.join(outputdir, cfile0)) and os.path.exists(gfile):
            nn1 = oapackage.nArrayFile(cfile)
            nn2 = oapackage.nArrayFile(gfile)
            raise Exception('both .oa and .oa.gz exist: %s' % cfile)
            if nn2 > nn1:
                print('### removing %s!!!' % cfile)
                os.remove(cfile)

        nn = oapackage.nArrays(cfile)
        mode = tagtype[ii]
        cfilex = oapackage.oahelper.checkOAfile(cfile)
        if cfilex is not None:
            cfilebase = os.path.basename(cfilex)
        else:
            cfilebase = None
        if nn >= 0:
            break

    if verbose:
        print('cdesignTag: N %d, kk %d: selected tag %s: nn %d' %
              (N, kk, tag, nn))
    # special case
    if kk == N and tag == 'cdesign-diagonal':
        mode = 'full'

    if verbose >= 2:
        print(cfile)
    return cfile, nn, mode


def generateConferenceResults(presults, ll, ct=None, full=None):
    pareto_results = presults
    pareto_results['type'] = 'conference designs'
    pareto_results['arrayfile'] = None
    pareto_results['presults'] = None

    pareto_results['full'] = full
    pareto_results['full_results'] = full
    pareto_results['idstr'] = 'cdesign-%d-%d' % (pareto_results['N'], pareto_results['ncolumns'])
    if ct is not None:
        pareto_results['ctidstr'] = ct.idstr()
        assert (ct.N == pareto_results['N'])
        assert (ct.ncols == pareto_results['ncolumns'])

    pareto_results['narrays'] = len(ll)

    pareto_results['pareto_designs'] = [np.array(array) for array in presults['pareto_designs']]
    return pareto_results


# %% Webpage generation


def nprevzero(N, k, ncache):
    """ Return true if any previous result was zero """
    for ix in range(k - 1, 2, -1):
        # print(ix)
        p = ncache['full'].get('N%dk%d' % (N, ix), -1)
        # p = ncache['full'].get('N%dk%d' % (N, ix), -1)
        if p == 0:
            return True
    return False


def htmlTag(nn, kk, N, mode='full', href=None, ncache=None, verbose=0):
    """ Create html tag for number of designs

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
        if mode == 'full':
            txt = '%d' % nn
        else:
            if nn == 0:
                txt = '?'
            else:
                txt = '&ge; %d' % nn
        if href is not None:
            if nn < 6000 and nn > 0 or (href.endswith('html')):
                ss = e.a(txt, href=href, style='text-decoration: none;')
                hyper_link = True
                txt = ss               
        else:
            pass
    else:
        if verbose:
            print('htmlTag: nn is negative')
        if kk <= N:
            if ncache is None:
                txt = '?'
            else:
                if verbose >= 1:
                    print('htmlTag: mode %s, N %d, k %d' % (mode, N, kk))
                if nprevzero(N, kk, ncache):
                    if verbose >= 1:
                        print('htmlTag: nprevzero(%d, %d, ncache) is True' % (N, kk))
                    txt = ''
                else:
                    txt = '?'
        else:
            txt = ''
    return txt, hyper_link


def latexResults(outputdir):
    X = []
    print('make latex results table...')

    NN = range(4, 31, 2)
    kk = range(0, np.max(NN) + 1)

    X = np.zeros((1 + 1 + len(kk) - 2, 1 + len(NN)), dtype=object)
    X[:] = ''
    X[0, 0] = ''
    X[1, 0] = '$k$'
    for ii, N in enumerate(NN):
        X[1, 1 + ii] = N
        X[0, 1 + ii] = ''
        for ki, k in enumerate(range(2, N + 1)):
            if k > N:
                X[1 + 1 + ki, 1 + ii] = ''
            else:
                cfile0 = 'cdesign-%d-%d.oa' % (N, k)
                nn = oapackage.nArrays(os.path.join(outputdir, cfile0))
                if nn < 0 and k == N:
                    cfile0 = 'cdesign-diagonal-%d-%d.oa' % (N, k)
                    nn = oapackage.nArrays(os.path.join(outputdir, cfile0))
                if nn < 0:
                    cfile0 = 'cdesign-diagonal-%d-%d.oa' % (N, k)
                    nnm = oapackage.nArrays(os.path.join(outputdir, cfile0))
                    if nnm > 0:
                        X[1 + 1 + ki, 1 + ii] = '$\ge %d$' % nnm
                    else:
                        X[1 + 1 + ki, 1 + ii] = '?'
                else:
                    X[1 + 1 + ki, 1 + ii] = nn
                X[1 + 1 + ki, 0] = '%d' % k
    X[0, 1] = '$N$'

    X[2, 0] = X[2, 0] + r'\rule{0pt}{2.9ex}'
    return X


def createConferenceDesignsPageHeader(page, makeheader, conference_class, ncolumns, full_results=False):
    xstr = 'C(%d, %d)' % (conference_class.N, ncolumns)
    xstrplain = xstr
    if makeheader:
        page.init(title="Class %s" % xstrplain,
                  css=('../oastyle.css'),
                  lang='en', htmlattrs=dict({'xmlns': 'http://www.w3.org/1999/xhtml', 'xml:lang': 'en'}),
                  header="<!-- Start of page -->",
                  htmlheader=oaresearch.research.oaCssStyle(addframe=True),
                  bodyattrs=dict({'style': 'padding-left: 3px;'}),
                  doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
                  metainfo=({'text/html': 'charset=utf-8', 'keywords': 'conference designs',
                             'robots': 'index, follow', 'description': 'conference designs'}),
                  footer="<!-- End of page -->")

    if full_results:
        page.h1('Conference designs %s ' % xstr)
    else:
        page.h1('Conference designs %s (<b>partial results</b>) ' % xstr)
    oap = e.a('Orthogonal Array package', href='../../software.html')
    pstr = 'This page contains information about conference designs. '
    pstr += 'The results have been generated with the %s.' % oap
    pstr += ' If you use these results, please cite the paper ' + \
            citation('conference', style='full') + '.'
    page.p(pstr)


def createConferenceDesignsPageResultsTable(page, pareto_results, verbose=0):
    full_results = pareto_results.get('full')
    if full_results:
        page.h2('Results')
    else:
        page.h2('Results (partial results)')

    page.table()
    page.tr(style='font-weight: bold;')
    page.td('Statistic', style='padding-right:30px;')
    page.td(('Results'), style='padding-right:8px;')
    page.tr.close()

    def simpleRow(a, b):
        page.tr(style='')
        page.td(a, style='padding-right:30px;')
        page.td(b, style='padding-right:8px;')
        page.tr.close()

    narrays = pareto_results['narrays']
    simpleRow('Number of non-isomorphic designs', str(pareto_results['narrays']))

    if narrays > 0:
        if 'ranksecondorder_min' in pareto_results:
            simpleRow('Minimum/Maximum rank of model matrix with 2FI and QE',
                      '%d/%d' % (pareto_results['ranksecondorder_min'], pareto_results['ranksecondorder_max']))
            simpleRow('Minimum/Maximum rank of model matrix for 2FI',
                      '%d/%d' % (pareto_results['rankinteraction_min'], pareto_results['rankinteraction_max']))
        else:
            simpleRow('Maximum rank of model matrix with 2FI and QE', '%d' % (pareto_results['ranksecondorder_max']))
            simpleRow('Maximum rank of model matrix for 2FI', '%d' % (pareto_results['rankinteraction_max']))
        if 'B4_min' in pareto_results:
            simpleRow('Minimum B4', '%.4f' % pareto_results['B4_min'])
        if 'B4_max' in pareto_results:
            simpleRow('Maximum B4', '%.4f' % pareto_results['B4_max'])
        if 'F4_min' in pareto_results:
            simpleRow('Minimum F4', '%s' % (pareto_results['F4_min'],))
        if 'F4_max' in pareto_results:
            simpleRow('Maximum F4', '%s' % (pareto_results['F4_max'],))
    if 'totaltime' in list(pareto_results.keys()):
        simpleRow('Processing time', str(pareto_results['totaltime']) + 's')

    if pareto_results.get('datafile_tag', None) is not None:
        simpleRow('Data', pareto_results['datafile_tag'])

    page.table.close()
    page.p(style='font-size: smaller;')
    page.add('Note: 2FI: two-factor interaction; QE: quadratic effect')
    page.p.close()


def createConferenceDesignsPageLoadDesignsFile(pareto_results, htmlsubdir=None, verbose=1):
    havearrayfile = 0
    if 'arrayfile' in list(pareto_results.keys()):
        if not pareto_results['arrayfile'] is None:
            havearrayfile = 1

    if verbose:
        print('createConferenceDesignsPageLoadDesignsFile: havearrayfile %d' % havearrayfile)
    if havearrayfile:
        iafile = pareto_results['arrayfile']
        outfile0 = pareto_results['idstr'] + '.oa'

        na = oapackage.nArrayFile(iafile)
        if verbose >= 2:
            print('conferenceDesignsPage: read arrayfile %s: na %d' % (iafile, na))

        if na < 5000 and na >= 0:
            if htmlsubdir is None:
                raise Exception('need html subdirectory to copy .oa file')
            outfilefinal = oaresearch.filetools.copyOAfile(iafile, htmlsubdir, outfile0, convert='T', zipfile=None,
                                                           verbose=1, cache=0)

            if pareto_results.get('full', False):
                if verbose:
                    print('conferenceDesignsPage: full results')
                pareto_results['datafilestr'] = 'all arrays'
                htag = oaresearch.htmltools.formatArrayHyperlink(
                    pareto_results['datafilestr'], outfile0, iafile)

                pareto_results['datafilestr'] = 'all arrays'
                pareto_results['datafile_tag'] = htag
            else:
                na = oapackage.nArrayFile(os.path.join(htmlsubdir, outfilefinal))
                pareto_results['datafilestr'] = '%d array(s)' % na
                htag = oaresearch.htmltools.formatArrayHyperlink(
                    pareto_results['datafilestr'], outfile0, iafile)
                pareto_results['datafile_tag'] = htag

        else:
            if verbose:
                print('conferenceDesignsPage: no datafile (na %d)' % na)
            pareto_results['datafilestr'] = '-'
            pareto_results['datafile_tag'] = None


def createConferenceDesignsPageParetoTable(page, pareto_results, verbose=0, htmlsubdir=None):
    """ Create table with Pareto results and add to the markup object
    
    Args:
        page (markup.page): html page to add table to
        
    Returns:
        rtable (array): generated table
    """
    if verbose:
        print('createConferenceDesignsPageParetoTable: start')

    pareto_indices = pareto_results['pareto_indices']
    pareto_data = pareto_results['pareto_data']

    if pareto_results['narrays'] > 0 and pareto_results.get('full_results'):
        add_extra = True
        if verbose:
            print('createConferenceDesignsPageParetoTable: do statistics2htmltable')

        header = ['Index Pareto file', 'Index design file', 'Rank 2FI and QE', 'Rank 2FI', 'F4', 'B4']
        if add_extra:
            for tag in ['PEC', 'PIC', 'PPC']:
                for kk in [4, 5]:
                    header += [tag + '%d' % kk]

        rtable = np.zeros((1 + len(pareto_results['pareto_indices']), len(header)), dtype='|U208')
        rtable[:] = ' '
        for ii, h in enumerate(header):
            rtable[0, ii] = header[ii]

        sort_indices = oapackage.sortrows(np.array([p['F4'] for p in pareto_results['pareto_data']]))

        for ii, sort_index in enumerate(sort_indices):
            pareto_idx = sort_index
            array_list_idx = pareto_indices[sort_index]
            rank_secondorder = str(pareto_data[pareto_idx]['ranksecondorder'])
            rank_interaction = str(pareto_data[pareto_idx]['rankinteraction'])
            rowdata = ['%d' % pareto_idx, '%d' % array_list_idx, rank_secondorder, rank_interaction,
                       str(pareto_data[pareto_idx]['F4']), '%.2f' % ((pareto_data[pareto_idx]['B4']))]
            rtable[ii + 1, 0:len(rowdata)] = rowdata
            column_offset = len(rowdata)
            if add_extra:
                for tag in ['PEC', 'PIC', 'PPC']:
                    for kk in [4, 5]:
                        rtable[ii + 1, column_offset] = '%.3f' % (pareto_data[pareto_idx][tag + '%d' % kk])
                        column_offset = column_offset + 1

        subpage = oaresearch.research.array2html(rtable, header=1, tablestyle='border-collapse: collapse;',
                                                 trclass='', tdstyle='padding-right:1em;', trstyle='',
                                                 thstyle='text-align:left; padding-right: 1em;', comment=None)
        page.br(clear='both')
        page.h2('Pareto optimal designs')
        page.p()
        if pareto_results['nclasses'] == 1:
            pareto_classes_text = 'in %d class' % pareto_results['nclasses']
        else:
            pareto_classes_text = 'in %d classes' % pareto_results['nclasses']
        if pareto_results['npareto'] == 1:
            page.add('There is %d Pareto optimal design %s.' % (pareto_results['npareto'], pareto_classes_text))
        else:
            page.add('There are %d Pareto optimal designs %s.' % (pareto_results['npareto'], pareto_classes_text))
        pareto_type = pareto_results['pareto_type']
        if ',' in pareto_type:
            k = pareto_type.rfind(", ")
            pareto_type = pareto_type[:k] + ", and " + pareto_type[k + 1:]

        page.add('Pareto optimality is according to %s (any other statistics are ignored).' % pareto_type)
        page.p.close()
        if pareto_results.get('pareto_designs', None) is not None:
            pdesigns = pareto_results.get('pareto_designs', None)

        pfile0 = pareto_results['idstr'] + '-pareto.oa'

        if htmlsubdir is not None:
            pfile = os.path.join(htmlsubdir, pfile0)
            oapackage.writearrayfile(pfile, [oapackage.array_link(array) for array in pdesigns])
            page.p('All %s' % e.a('Pareto optimal designs', href=pfile0) + '.')

        page.add(str(subpage))
        page.br(clear='both')
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
    latextable = oapackage.array2latex(rtable_latex, hlines=[0],
                                       comment=['conference desgins N=%d, k=%d' % (N, ncolumns),
                                                'offset for indices is 1'])
    return latextable


def conferenceDesignsPage(pareto_results, verbose=1, makeheader=True, htmlsubdir=None, generate_latex=True, html_template=False):
    """ Generate html page for class conference designs
    
    Args:
        pareto_results (dict): structure with results for Pareto optimal designs
        html_template (bool): If True then place html files in subfolder templates/
    """

    N = pareto_results['N']
    ncolumns = pareto_results['ncolumns']

    if makeheader:
        pext = 'html'
    else:
        pext = 'html.template'

    conference_class = oapackage.conference_t(N, ncolumns, 0)

    if verbose:
        print('conferenceDesignsPage: start of generation')

    page = markup.page()
    createConferenceDesignsPageHeader(page, makeheader=makeheader, conference_class=conference_class, ncolumns=ncolumns,
                                      full_results=pareto_results['full_results'])
    createConferenceDesignsPageLoadDesignsFile(pareto_results, htmlsubdir=htmlsubdir)
    createConferenceDesignsPageResultsTable(page, pareto_results, verbose=verbose)
    rtable = createConferenceDesignsPageParetoTable(page, pareto_results, verbose=verbose, htmlsubdir=htmlsubdir)

    latextable = _convert_to_latex_table(rtable, N, ncolumns)

    if generate_latex:
        with open(os.path.join(htmlsubdir, 'conference-N%dk%d.tex' % (N, ncolumns)), 'wt') as fid:
            fid.write(latextable)
        import pickle
        with open(os.path.join(htmlsubdir, 'conference-N%dk%d-rtable.pickle' % (N, ncolumns)), 'wb') as fid:
            pickle.dump({'rtable': rtable, 'N': pareto_results['N'], 'ncolumns': pareto_results['ncolumns']}, fid)

    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)

    page.p('<br/>\n')
    page.p('Page generated on %s.' % dstr)

    pstr = str(page).replace('<meta content="charset=utf-8" name="text/html" />',
                             '<meta http-equiv="Content-Type" content="charset=utf-8" name="text/html" />')

    return pstr
