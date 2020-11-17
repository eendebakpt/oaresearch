# -*- coding: utf-8 -*-
"""

Example script for calculating conference designs

Pieter Eendebak <pieter.eendebak@gmail.com>
"""


# %% Load necessary packages
from typing import List, Any
import tempfile
import os
import sys
import platform
import numpy as np
import time
from importlib import reload
from os.path import join
import pickle
import pdb
import shutil
import webbrowser

import oapackage
from oapackage.oahelper import write_text_arrayfile
import oapackage.graphtools
from oapackage.markup import oneliner as e
from oaresearch.research_conference import htmlTag
import oaresearch.research_conference
from oaresearch.research_conference import calculateConferencePareto, conferenceResultsFile, generateConferenceResults, \
    conferenceDesignsPage
from oapackage.conference import conferenceProjectionStatistics
from oaresearch.research_conference import SingleConferenceParetoCombiner, generate_conference_latex_tables
from oapackage import markup
from oapackage.oahelper import create_pareto_element
from oaresearch.research_conference import generate_or_load_conference_results, createConferenceParetoElement, \
    calculateConferencePareto, createConferenceDesignsPageParetoTable, cdesignTag, conference_design_has_extensions
from oaresearch.research_conference import select_even_odd_conference_designs, generate_even_odd_conference_designs

generate_webpage = True


# %% Setup directories
resultsdir = join(os.path.expanduser('~'), 'oatmp')

# outputdir of the code generating the designs
outputdir = oapackage.mkdirc(
    os.path.join(os.path.expanduser('~'), 'oatmp', 'conf'))

if platform.node() == 'woelmuis':
    paperdir = '/home/eendebakpt/misc/oa/article-conference/'

if generate_webpage:
    htmldir = os.path.join(os.path.expanduser('~'), 'misc', 'oapage2')
    html_template = True
    if 0:
        # for testing...
        htmldir = os.path.join(os.path.expanduser('~'), 'oatmp', 'confpage_dc')
        html_template = False
    oapackage.mkdirc(os.path.join(htmldir))
    conference_html_dir = oapackage.mkdirc(os.path.join(htmldir, 'conference'))


def conference_designs_result_file(outputdir: str, tag: str, N: int, kk: int) -> str:
    cache_tag = f'result-pareto-v1-{oaresearch.research_conference.conferenceParetoIdentifier()}'
    oapackage.mkdirc(os.path.join(outputdir, cache_tag))
    cachefile = os.path.join(
        outputdir, cache_tag, tag + '-' + 'N%dk%d' % (N, kk) + '.pickle')
    return cachefile


# %%
if 0:
    addMaximumExtensionColumns = True
    N = 24
    kk = 16
    N = 22
    kk = 8
    # N=20; kk=13;
    N = 12
    kk = 4
    #N = 28;kk = 8
    # N = 4; kk = 2
    # N=30;kk=6

    t0 = time.time()

    cfile, nn, mode = conferenceResultsFile(N, kk, outputdir, tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
                                            tagtype=['full', 'r', 'r'], verbose=1)
    ll = oapackage.readarrayfile(cfile)
    ll = ll[0:]

    pareto_results, cfile = generate_or_load_conference_results(N, kk, outputdir, dc_outputdir=resultsdir,
                                                                double_conference_cases=[10, 16, 24], addMaximumExtensionColumns=addMaximumExtensionColumns)
    page = conferenceDesignsPage(pareto_results, verbose=1, makeheader=True,
                                 htmlsubdir=conference_html_dir, html_template=html_template)
    dt = time.time() - t0
    print('processing time: %.1f [s]' % dt)

    oapackage.oahelper.testHtml(str(page))

    # 600 seconds for N=20, kk=13
    # with refactoring and mkl: 251 [s]

    # N24k17: 3 sec/array full,

    # TODO: generate data for C(24, k)
    # TODO: generate page for C(24, k)

    for ii, al in enumerate(ll):
        print('design %d' % ii)
        sys.stdout.flush()
        conference_design_has_extensions(al, verbose=1)

# %%
if 0:
    reload(oaresearch.research_conference)
    from oaresearch.research_conference import SingleConferenceParetoCombiner

    cache_dir = oapackage.mkdirc(os.path.join(
        resultsdir, 'doubleconference-%d' % (2 * N), 'sc_pareto_cache'))
    pareto_calculator = SingleConferenceParetoCombiner(
        outputdir, cache_dir=cache_dir, cache=True)
    presults = pareto_calculator.load_combined_results(5)
    print(presults)

# %%

# %%
if 0:
    # designs = [oapackage.array_link(al) for al in pareto_results['pareto_designs']]

    # designs=ll[5168:5171]
    # designs=ll[4347:4349]

    al1 = ll[4347]
    al2 = ll[4506]

    designs = ll[4000:]
    for jj, al in enumerate(designs):
        f4, b4, rank, rankq = oaresearch.research_conference.conferenceStatistics(
            al, verbose=0)
        pec5, pic5, ppc5 = conferenceProjectionStatistics(al, 5)
        if np.abs(ppc5 - 1.2130420253) < 1e-4:
            pec, pic, ppc = conferenceProjectionStatistics(al, 4)
            print('array %d: b4 %s f4 %s' % (jj, b4, f4))
            print('  pec %s, pec5 %s, ppc %s, ppc5 %s' %
                  (pec, pec5, ppc, ppc5))

# %%
if 0:
    rtable = createConferenceDesignsPageParetoTable(
        markup.page(), pareto_results, verbose=2, htmlsubdir=None)
    latextable = oapackage.array2latex(rtable)
    print(latextable)

# %%

if 0:
    from oapackage.oahelper import create_pareto_element
    from oaresearch.research_conference import createConferenceParetoElement, calculateConferencePareto

    presults = calculateConferencePareto(ll, N=None, k=None, verbose=1)


# %% Generate subpages for the designs


def conferenceSubPages(tag='conference', Nmax=40, Nstart=4, kmax=None, outputdir=None, conference_html_dir=None,
                       verbose=1, specials={}, Nstep=2, NmaxPareto=40, cache=True,
                       double_conference_cases=(24,), html_template=False, addMaximumExtensionColumns=False):
    """ Generate subpages for single conference results

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
        print('conferenceSubPages: tag %s, Nmax %d, kmax %d' %
              (tag, Nmax, kmax))

    subpages = {}
    subpages[tag] = {}

    for ki, kk in enumerate(krange):
        for Ni, N in enumerate(Nrange):

            if tag == 'conference':
                if kk > N:
                    continue

            subpages[tag]['N%dk%d' % (N, kk)] = {}
            cachefile = conference_designs_result_file(outputdir, tag, N, kk)

            if cache and os.path.exists(cachefile):
                with open(cachefile, 'rb') as fid:
                    if verbose >= 2:
                        print('loading results from cachefile %s' % cachefile)
                    pareto_results, cfile = pickle.load(fid)
                if verbose:
                    print(
                        'conferenceSubPages %s: from cache: N %d, columns %d' % (tag, N, kk))
                if pareto_results['presults'] is not None:
                    pversion = pareto_results['_version']
                    cversion = oaresearch.research_conference.conferenceParetoIdentifier()
                    if pareto_results['_version'] != cversion:
                        raise Exception(
                            f'conference Pareto definition was updated! file {cachefile} data {pversion} code {cversion}')

            else:
                pareto_results, cfile = generate_or_load_conference_results(N, kk, outputdir, dc_outputdir=resultsdir,
                                                                            double_conference_cases=double_conference_cases, addMaximumExtensionColumns=N <= 20)

                if verbose >= 2:
                    print('storing results in cachefile %s' % cachefile)

                with open(cachefile, 'wb') as fid:
                    pareto_results['presults'] = None
                    pickle.dump((pareto_results, cfile), fid)

            if pareto_results.get('narrays', None) is None:
                continue

            # create HTML page
            page = conferenceDesignsPage(pareto_results, verbose=verbose >= 2, makeheader=True,
                                         htmlsubdir=conference_html_dir, html_template=html_template)

            # write results
            htmlfile0 = tag + '-%d-%d.html' % (N, kk)
            if html_template:
                oapackage.mkdirc(os.path.join(
                    conference_html_dir, 'templates'))
                htmlfile = os.path.join(
                    conference_html_dir, 'templates', htmlfile0 + '.template')
            else:
                htmlfile = os.path.join(conference_html_dir, htmlfile0)

            sx = subpages[tag]['N%dk%d' % (N, kk)]
            sx['htmlpage0'] = htmlfile0
            sx['htmlpage'] = htmlfile
            sx['pareto_results'] = pareto_results
            sx['arrayfile'] = cfile

            if verbose >= 2:
                print('writing to %s' % htmlfile)
            with open(htmlfile, 'wt') as fid:
                fid.write(page)

    conferenceSubPages._debugdata = debugdata
    return subpages


# generated_subpages = conferenceSubPages(tag='cdesign', Nmax=40, Nstart=4, verbose=2, cache=False, outputdir=outputdir)
generated_subpages = conferenceSubPages(tag='cdesign', Nmax=40, Nstart=4, verbose=1, cache=True,
                                        outputdir=outputdir, conference_html_dir=conference_html_dir, html_template=html_template)
#generated_subpages = conferenceSubPages(tag='cdesign', Nmax=24, Nstart=4, verbose=1, kmax=9, cache=True, outputdir=outputdir, conference_html_dir = conference_html_dir, html_template = html_template)

# %% Results table for latex

htmlsubdir = os.path.join(htmldir, 'conference')

generate_conference_latex_tables(htmlsubdir, verbose=1)

#%%

#%% Determine number of extensible designs
from oaresearch.research_conference import conference_design_extensions


for Nx in range(4, 10, 2):
  for k in range(2, Nx+1):
    data=generated_subpages['cdesign'][f'N{Nx}k{k}' ]
    
    data.keys()
    
    pr=data['pareto_results']
    data['arrayfile']
    
    cfile=os.path.join(outputdir, pr['idstr']+'.oa' )
    #print(cfile)
    
    designs=oapackage.readarrayfile(cfile)
    
    ee=[len(conference_design_extensions(design))>0 for design in designs]
    #print(ee)
    ne = np.sum(ee)
    print(f'N {Nx} k {k}: {len(designs)} designs, {ne} can be extended')


# %% Testing
if 1:
    # subpage='xxx.html'

    tag = 'cdesign'
    page = markup.page()
    N = 12
    kk = 9
    generated_result = generated_subpages['cdesign'].get(
        'N%dk%d' % (N, kk), None)
    subpage = generated_result['htmlpage0']
    cdesignTag(N, kk, page, outputdir, '',
               ['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
               ['full', 'r', 'r'], verbose=2, subpage=subpage, generated_result=generated_result, conference_html_dir=conference_html_dir)

    generated_result = generated_subpages['cdesign'].get(
        'N%dk%d' % (28, 4), None)
    cdesignTag(28, 4, page, outputdir, '',
               ['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
               ['full', 'r', 'r'], verbose=2, subpage=subpage, generated_result=generated_result, conference_html_dir=conference_html_dir)

    tag = 'dconferencej1j3'
    N = 80
    kk = 8

    cdesignTag(N, kk, page, outputdir, tdstyle='', tags=[tag, tag + '-r'],
               tagtype=['full', 'r'], verbose=2, ncache=None, subpage=None, generated_result=None, conference_html_dir=conference_html_dir)
    # tag = 'dconferencej1'
    # cdesignTag(N=N, kk=kk, page=page, outputdir=outputdir,
    # tags=[tag, tag + '-r'], tagtype=['full', 'r'], verbose=2, ncache=ncache)


# %%


def specialData():
    """ Special data """
    s = dict()
    s['cdesign'] = {}

    # generated by Alan
    tt = {8: 1839474,
          9: 8259167,
          10: 8667156,
          11: 4124471,
          12: 2397144,
          13: 1806230,
          14: 1353790,
          15: 888475,
          16: 499614,
          17: 234006,
          18: 91773,
          19: 28730,
          20: 7417,
          21: 1377,
          22: 232,
          23: 19,
          # 24: 9,
          }
    for k in tt:
        tt[k] = (tt[k], 'full')
    # s['cdesign'][24] = tt

    # generated by Alan
    tt = {  # 5: 14,
        6: 209,
        7: 897,
        8: 6769,
        9: 48980,
        10: 107503,
        11: 65517,
        12: 16184,
        13: 2633,
        14: 1088,
        15: 675,
        16: 455,
        17: 295,
        18: 221,
        19: 168,
        20: 132,
        21: 92,
        22: 76,
        23: 63,
        24: 57,
        25: 50,
        26: 48,
        27: 44,
        28: 41}
    for k in tt:
        tt[k] = (tt[k], 'r')
        if k == 28:
            tt[k] = (tt[k][0], 'full')

    s['cdesign'][28] = tt

    # from doubleconference_split_analyse
    ttx = {5: [30, 8],
           6: [1588, 103],
           7: [87929, 762],
           8: [1839474, 1215],
           9: [8259167, 110],
           10: [8667156, 0],
           11: [4124471, 0],
           12: [2397144, 0],
           13: [1806230, 0],
           14: [1353790, 0],
           15: [888475, 0],
           16: [499614, 0],
           17: [234006, 0],
           18: [91773, 0],
           19: [28730, 0],
           20: [7417, 0],
           21: [1377, 0],
           22: [232, 0],
           23: [19, 0],
           24: [9, 0]}
    tt = {}
    tteo = {}
    for k in ttx:
        if k <= 7:
            continue
        tt[k] = (ttx[k][0] + ttx[k][1], 'full')
        tteo[k] = (ttx[k][1], 'full')
    s['dconferencej1j3'] = {48: tt}
    s['dconferencej1j3-eo'] = {48: tteo}
    return s


specialdata = specialData()
# specialdata = {}

specialdataDC = specialdata


# %%


def DconferencePage(page, tag='dconference', Nmax=26, Nstart=4, kmax=None,
                    ta='left', verbose=1, specials={}, Nstep=2,
                    tableclass='conftable', tdstyle=None, subpages=None, conference_html_dir=None):
    """ Generate a table with matrices

    Arguments:
        page (markup object)
        tag (str): type of data
        Nmax, Nstart, kmax (int)
        ta (str): text-align property
        verbose (int): verbosity level

    """

    if tdstyle is None:
        if ta == 'left':
            tar = 'right'
        else:
            tar = 'left'
        tdstyle = 'text-align:%s; margin-%s: 1em; padding-%s:1em; margin-%s: 1px;' % (
            ta, tar, tar, ta)

    Nrange = range(Nstart, Nmax + 1, Nstep)
    if kmax == -2:
        kmax = int(np.ceil(Nmax / 2) + 1)
    if kmax is None:
        kmax = np.max(Nrange) + 2
    krange = range(2, kmax)
    ncols = len(Nrange) + 1
    if verbose:
        print('Dconferencepage: tag %s, Nmax %d, kmax %d' % (tag, Nmax, kmax))
    page.p()
    page.table(class_=tableclass)
    for ii in range(ncols):
        page.col()
    page.tr()
    page.th('', style='color: white;')
    page.th(['Number of rows'], colspan=len(
        Nrange) + 0, style=tdstyle)
    page.tr.close()
    page.tr()
    page.th(['k'] + ['%d' % k for k in Nrange], style=tdstyle)
    page.tr.close()

    ncache = {}
    for ki, kk in enumerate(krange):
        page.tr()
        page.td('%s' % str(kk), style=tdstyle)
        for Ni, N in enumerate(Nrange):
            subpage = None
            generated_result = None
            if subpages is not None:
                generated_result = subpages.get('N%dk%d' % (N, kk), None)
                if generated_result is not None and len(generated_result) != 0:
                    subpage = generated_result['htmlpage0']

            if tag in specials:
                if N in specials[tag]:
                    if kk in specials[tag][N]:

                        print('special case: %s N %d, kk %d' % (tag, N, kk))
                        nn, mode = specials[tag][N][kk]
                        if 1:
                            txt, link = htmlTag(nn, kk, N, mode=mode)
                            page.td(txt, style=tdstyle)

                            continue
            if verbose:
                print('DconferencePage: tag %s: %d %d: subpage %s' %
                      (tag, N, kk, subpage))

            if tag == 'cdesign':
                cdesignTag(N, kk, page, outputdir, tdstyle,
                           ['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
                           ['full', 'r', 'r'], subpage=subpage, generated_result=generated_result, conference_html_dir=conference_html_dir)
            else:
                cdesignTag(N, kk, page, outputdir, tdstyle,
                           tags=[tag, tag + '-r'], tagtype=['full', 'r'], ncache=ncache, conference_html_dir=conference_html_dir)

        page.tr.close()
    page.table.close()
    page.p.close()

    return {'ncache': ncache}


# %% Get even-odd designs


if 0:
    generate_even_odd_conference_designs(outputdir=outputdir)


# %%


if generate_webpage:
    tdstyle = 'text-align:left; padding-left: .6em; margin-right: 0em; padding-right: .6em; margin-left: 0px;'

    citation = oaresearch.research.citation('conference', style='brief')
    citation_enumeration = oaresearch.research.citation(
        'cenumeration', style='brief')
    page = oapackage.markup.page()

    if not html_template:
        page.init(title="Conference designs",
                  css=('oastyle.css'),
                  lang='en',
                  header="", htmlattrs=dict({'xmlns': 'http://www.w3.org/1999/xhtml'}),
                  # doctype=markup.doctype.strict,
                  doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"   "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">',
                  metainfo=({'keywords': 'conference designs', 'robots': 'index, follow',
                             'description': 'Complete Enumeration of conference designs'}),
                  footer="<!-- End of page -->")

    page.h1('Non-isomorphic conference designs', style='margin-bottom: 0.1em;')
    page.p('<em style="color:darkblue;">pieter.eendebak@gmail.com</em>',
           style='margin-top: 0px;')
    ss = 'On this page we present numbers of isomorphism classes for conference designs with a specified number of runs (or rows) and factors (or columns). '
    ss += 'For all the cases a set of representatives for the isomorphism classes is available. The algorithm used to generate the results is described in %s and %s. ' % (
        citation, citation_enumeration)
    page.p(ss)

    if 1:
        page.h2('Single conference designs')
        page.p()
        page.add(
            'Conference designs are matrixes of size 2m x k with values 0, +1 or -1, where m &ge; 2 and k &le; m. Each column contains exactly one zero and each row conains at most one zero. The columns are orthogonal to each other.')
        ss = 'Square single conference designs are <a href="https://en.wikipedia.org/w/index.php?title=Conference_matrix&oldid=702633334">conference matrices</a>.'
        page.add(ss)
        page.p.close()

        page.p(
            'A definitive screening design can be constructed by appending a conference design with both its fold-over and a row of zeroes.')
        DconferencePage(page, tag='cdesign', Nmax=41,
                        kmax=32, specials=specialdata, tdstyle=tdstyle,
                        tableclass='conftable conftable1', subpages=generated_subpages['cdesign'], conference_html_dir=conference_html_dir)

    if 1:
        page.h2('Double conference designs (DCDs)')
        page.p()
        page.add(
            'Double conference designs are matrices of size 4m x k with values 0, +1 or -1, where m &ge;2 and k &le; 2m. Each column contains exactly two zeros and each row contains at most one zero. The columns are orthogonal to each other.')
        page.p.close()

        def subheader(s):
            """ Create subheader in html page """
            page.h3(s)

        subheader('DCDs with level balance and orthogonal interaction columns.')
        page.p(
            'In double conference designs with level balance and orthogonal interaction columns, the elements of each column sums to zero. In addition, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero too. A definitive screening design can be constructed by appending a double conference design with a row of zeroes.')
        # page.p('All J1 and J3 values are zero.')
        DconferencePage(page, tag='dconferencej1j3', specials=specialdata,
                        Nstart=4, Nmax=81, kmax=27, Nstep=4, tableclass='conftable conftable2', conference_html_dir=conference_html_dir)
        page.p()

        if 0:
            page.add('Only the even-odd designs')
            page.p.close()
            DconferencePage(page, tag='dconferencej1j3-eo',
                            Nstart=4, Nmax=81, kmax=21, Nstep=4, specials=specialdata, tableclass='conftable conftable2', conference_html_dir=conference_html_dir)

        if 1:
            subheader('DCDs with level balance')  # J1 values are zero
            page.p(
                'In double conference designs with level balance, the elements of each column sum to zero.')
            r = DconferencePage(page, tag='dconferencej1', kmax=14, tableclass='conftable conftable3',
                                conference_html_dir=conference_html_dir)

        if 1:
            subheader('DCDs with orthogonality of interaction columns')  # J3=0
            page.p(
                'In double conference designs with orthogonal interaction columns, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero.')
            DconferencePage(page, tag='dconferencej3', kmax=14, Nmax=24,
                            tableclass='conftable conftable4', conference_html_dir=conference_html_dir)

        if 1:
            subheader('Plain DCDs')  # no restrictions on J1 and J3
            page.p(
                'In double conference designs with orthogonal interaction columns, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero.')
            DconferencePage(page, tag='dconference2', Nmax=20, kmax=13,
                            tableclass='conftable conftable5', conference_html_dir=conference_html_dir)

        page.h2('Weighing matrices')

        if 1:
            # no restrictions on J1 and J3
            subheader('Plain DCDs (full isomorphism class)')
            page.p()
            page.add(
                'For the isomorphism class level-permutatations of the rows are allowed.')
            page.add(
                'We also allow multiple zeros in a row.')
            page.add(
                'The square double conference designs of this type form a complete non-isomorphic set of weighing matrices of type W(N, N-2).')
            page.add(
                'The square double conference matrices are a complete non-isomorphic set of weighing matrices of type W(N, N-2).')
            page.p.close()
            DconferencePage(page, tag='weighing2', kmax=15,
                            Nmax=25, conference_html_dir=conference_html_dir)

    ss = 'If you make use of these results, please cite the paper %s.' % citation
    page.p(ss)
    if 1:
        ss = 'The square conference designs are <a href="https://en.wikipedia.org/w/index.php?title=Conference_matrix">conference matrices</a>.'
        page.p(ss)

    if 0:
        page.h1('Obsolete results')
        subheader('DCDs with level balance and orthogonal interaction columns.')
        page.p(
            'These designs have N &equiv; 2 mod 4. In double conference designs with level balance and orthogonal interaction columns, the elements of each column sums to zero. In addition, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero too. A definitive screening design can be constructed by appending a double conference design with a row of zeroes.')
        DconferencePage(page, tag='dconferencej1j3', Nstart=6, Nmax=81, kmax=28,
                        Nstep=4, conference_html_dir=conference_html_dir)

        subheader('Double conference matrices (no restrictions on J1 and J3)')
        page.p()
        page.add(
            'Double conference matrices are matrixes of size 2m x k. Each columns contains exactly two zeros and an arbitrary number of values +1 or -1. The columns are orthogonal to each other.')
        page.add(
            'The square double conference matrices are weighing matrices of type W(N, N-2).')
        page.p.close()
        DconferencePage(page, tag='dconference2',
                        conference_html_dir=conference_html_dir)

    if html_template:
        hfile = os.path.join(htmldir, 'templates', 'conference.html')
    else:
        hfile = os.path.join(htmldir, 'conference.html')

    with open(hfile, 'w') as fid:
        _ = fid.write(str(page))
    webbrowser.open_new_tab(hfile)


# %%


# designs=[oapackage.exampleArray(3)]
# filename=tempfile.mktemp(suffix='.oa')
#write_text_arrayfile(filename, designs, comment='Test')
# oapackage.readarrayfile(filename)

# %%
def createConferenceDesignsPageParetoTableDSD(
        page, pareto_results, verbose=0, htmlsubdir=None):
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

    add_extension_information = False
    add_maximum_extension_size = False

    if len(pareto_data) > 0:
        if pareto_data[0].get('has_extensions', None) is not None:
            add_extension_information = True
        if pareto_data[0].get('maximum_extension_size', None) is not None:
            add_maximum_extension_size = True

    ncolumns = pareto_results['ncolumns']

    if pareto_results['narrays'] > 0 and pareto_results.get('full_results'):
        add_extra = True

        header = ['k', 'Index design file',
                  'r_2FIandQE', 'r<sub>2FI</sub>', 'F<sub>4</sub>', 'B<sub>4</sub>']
        if add_extra:
            for tag in ['PEC', 'PIC', 'PPC']:
                for kk in [4, 5]:
                    header += [tag + '<sub>%d</sub>' % kk]
        if add_extension_information:
            header += ['Extensions']
        if add_maximum_extension_size:
            header += ['Max. columns']

        rtable = np.zeros(
            (1 + len(pareto_results['pareto_indices']), len(header)), dtype='|U208')
        rtable[:] = ' '
        for ii, h in enumerate(header):
            rtable[0, ii] = header[ii]

        sort_indices = oapackage.sortrows(
            np.array([p['F4'] for p in pareto_results['pareto_data']]))

        for ii, sort_index in enumerate(sort_indices):
            pareto_idx = sort_index
            array_list_idx = pareto_indices[sort_index]
            rank_secondorder = str(pareto_data[pareto_idx]['ranksecondorder'])
            rank_interaction = str(pareto_data[pareto_idx]['rankinteraction'])
            rowdata = [f'{ncolumns}', '%d' % array_list_idx, rank_secondorder, rank_interaction,
                       str(pareto_data[pareto_idx]['F4']), '%.2f' % ((pareto_data[pareto_idx]['B4']))]
            rtable[ii + 1, 0:len(rowdata)] = rowdata
            column_offset = len(rowdata)
            if add_extra:
                for tag in ['PEC', 'PIC', 'PPC']:
                    for kk in [4, 5]:
                        rtable[ii + 1, column_offset] = '%.3f' % (
                            pareto_data[pareto_idx][tag + '%d' % kk])
                        column_offset = column_offset + 1
            if add_extension_information:
                rtable[ii + 1, column_offset] = 'Yes' if pareto_data[pareto_idx]['has_extensions'] > 0 else 'No'
                column_offset = column_offset + 1
            if add_maximum_extension_size:
                rtable[ii + 1, column_offset] = pareto_data[pareto_idx]['maximum_extension_size']
                column_offset = column_offset + 1

    else:
        rtable = np.zeros((0, 0))

    return rtable


rtable = createConferenceDesignsPageParetoTable(None, pareto_results)

# %% Generate Pareto table
pdir = tempfile.mkdtemp(prefix=f'D_X_k')
pdir = r'/home/eendebakpt/Dropbox/conference designs/tmp/pareto-overview'

tag = 'cdesign'
citation_enumeration = oaresearch.research.citation(
    'cenumeration', style='brief')
page = markup.page()


page.init(title=f"Pareto optimal designs for D[*](N, k)",
          css=('oastyle.css'),
          lang='en',
          header="", htmlattrs=dict({'xmlns': 'http://www.w3.org/1999/xhtml'}),
          doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"   "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">',
          metainfo=({'keywords': 'conference designs', 'robots': 'index, follow',
                     'description': 'Complete Enumeration of conference designs'}),
          footer="<!-- End of page -->")

page.h1(f"Pareto optimal designs for D[*](N, k)", style='margin-bottom: 0.1em;')
page.p('<em style="color:darkblue;">pieter.eendebak@gmail.com</em>',
       style='margin-top: 0px;')
ss = f'On this page we present properties of the Pareto optimal designs in D(N, k). The optimality criteria are *,*,* and are discussed in detail in {citation_enumeration}. '
page.p(ss)

for N in [8, 10, 12, 14, 16, 18, 20]:
    DN = 2*N+1

    page.h3(f"Pareto optimal designs for D[*]({DN}, k)", style='margin-bottom: 0.1em;')

    rtables = []
    designs = []
    for number_of_columns in range(4, N+1):
        print(f'N {N}, properties of {number_of_columns} columns')

        cachefile = conference_designs_result_file(outputdir, tag, N, number_of_columns)

        with open(cachefile, 'rb') as fid:
            pareto_results, cfile = pickle.load(fid)

        subtable = createConferenceDesignsPageParetoTableDSD(None, pareto_results)
        rtables.append(subtable[1:, :])
        header = subtable[:1, :]
        designs += pareto_results['pareto_designs']

    designs = [oapackage.makearraylink(d) for d in designs]
    dsd_designs = [oapackage.conference2DSD(d) for d in designs]
    designfile_base = os.path.join(f'pareto-designs-N{N}.oa')
    dfile = os.path.join(pdir, designfile_base)
    write_text_arrayfile(dfile, dsd_designs, comment=f' Pareto optimal DSDs for {N} rows\nSee {citation_enumeration}')
    oapackage.readarrayfile(dfile)

    rtable = np.vstack(rtables)

    # for ii in [4,5]:
    #    header[0,2+ii]=f'PEC<sub>{ii}</sub>'

    rtable = np.vstack((header, rtable))
    rtable = np.hstack((rtable[:, :1], rtable[:, 2:]))

    subpage = oaresearch.research.array2html(rtable, header=1, tablestyle='border-collapse: collapse;',
                                             trclass='', tdstyle='padding-right:1em;', trstyle='',
                                             thstyle='text-align:left; padding-right: 1em;', comment=None)

    page.add(str(subpage))

    link = markup.oneliner.a(f'All Pareto DSDs with {DN} rows', href=designfile_base)
    page.p(f'{link} (ordered according the table above)')

page.p('TODO: check criteria with the ones used in the paper')

# %
hfile = os.path.join(pdir, 'pareto-page.html')
with open(hfile, 'w') as fid:
    _ = fid.write(str(page))
webbrowser.open_new_tab(hfile)

# %%
