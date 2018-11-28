# -*- coding: utf-8 -*-
"""

Example script for calculating conference designs

Pieter Eendebak <pieter.eendebak@gmail.com>
"""

#%% Load necessary packages
import os
import sys
import platform
import numpy as np
import time
from imp import reload
from os.path import join
import pickle
import pdb
import shutil

import oapackage
import oapackage.graphtools
from oapackage.markup import oneliner as e
from oaresearch.research_conference import htmlTag, nprevzero
import oaresearch.research_conference
reload(oaresearch.research_conference)
from oaresearch.research_conference import calculateConferencePareto, conferenceResultsFile, generateConferenceResults, conferenceDesignsPage
from oapackage.conference import conferenceProjectionStatistics
from oaresearch.research_conference import SingleConferenceParetoCombiner

generate_webpage = True

#%% Setup directories
resultsdir = join(os.path.expanduser('~'), 'oatmp')

# outputdir of the code generating the designs
outputdir = oapackage.mkdirc(
    os.path.join(os.path.expanduser('~'), 'oatmp', 'conf'))


if platform.node() == 'woelmuis':
    paperdir = '/home/eendebakpt/misc/oa/article-conference/'

if generate_webpage:
    htmldir = os.path.join(os.path.expanduser('~'),  'misc', 'oapage2')
    htemplate = True
    if 1:
        # for testing...
        htmldir = os.path.join(os.path.expanduser('~'),  'oatmp', 'confpage')
        htemplate = False
    oapackage.mkdirc(os.path.join(htmldir))
    cdir = oapackage.mkdirc(os.path.join(htmldir, 'conference'))


#%%
from oapackage import markup
from oapackage.oahelper import create_pareto_element
reload(oaresearch.research_conference)
from oaresearch.research_conference import createConferenceParetoElement, calculateConferencePareto, generateConferenceResults, conferenceDesignsPage, createConferenceDesignsPageParetoTable

N = 24; kk = 16
N = 12; kk = 4
#N=20; kk=13;
#N=12; kk=4;
N=16;kk=5
t0 = time.time()
cfile, nn, mode = conferenceResultsFile(N, kk, outputdir, tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'], tagtype=['full', 'r', 'r'], verbose=1)

ll = oapackage.readarrayfile(cfile)
ll = ll[0:]

from_double_conference = 1
if from_double_conference:
    dc_dir = os.path.join(resultsdir, 'doubleconference-%d' % (2*N))
    cache_dir = oapackage.mkdirc(os.path.join(dc_dir, 'sc_pareto_cache'))
    pareto_calculator = SingleConferenceParetoCombiner(dc_dir, cache_dir=cache_dir, cache=True)
    pareto_results=pareto_calculator.load_combined_results(kk)
    pareto_results['narrays']=pareto_results['combined_statistics']['number_conference_arrays']
    pareto_results['idstr'] = 'cdesign-%d-%d' % (pareto_results['N'], pareto_results['ncolumns'])

    pareto_results['full']=True
    pareto_results['full_results']=True
else:
    presults, pareto = calculateConferencePareto(ll, N=N, k=kk, verbose=1, addProjectionStatistics=None)
    pareto_results = generateConferenceResults(presults, ll, ct=None, full=mode == 'full')
    pareto_results['arrayfile'] = cfile

page = conferenceDesignsPage(pareto_results, verbose=1, makeheader=True, htmlsubdir=cdir)
dt = time.time() - t0
print('processing time: %.1f [s]' % dt)

oapackage.oahelper.testHtml(str(page))

# 600 seconds for N=20, kk=13
# with refactoring and mkl: 251 [s]

# N24k17: 3 sec/array full, 


# TODO: generate data for C(24, k)
# TODO: generate page for C(24, k)

#%%
reload(oaresearch.research_conference)
from oaresearch.research_conference import SingleConferenceParetoCombiner

cache_dir = oapackage.mkdirc(os.path.join(resultsdir, 'doubleconference-%d' % (2*N), 'sc_pareto_cache'))
pareto_calculator = SingleConferenceParetoCombiner(outputdir, cache_dir=cache_dir, cache=True)
presults=pareto_calculator.load_combined_results(5)
print(presults)

#%%

dc48 = os.path.join('resultsdir', '')
#%%
if 0:
    #designs = [oapackage.array_link(al) for al in pareto_results['pareto_designs']]

    # designs=ll[5168:5171]
    # designs=ll[4347:4349]

    al1 = ll[4347]
    al2 = ll[4506]

    designs = ll[4000:]
    for jj, al in enumerate(designs):
        f4, b4, rank, rankq = oaresearch.research_conference.conferenceStatistics(al, verbose=0)
        pec5, pic5, ppc5 = conferenceProjectionStatistics(al, 5)
        if np.abs(ppc5 - 1.2130420253) < 1e-4:
            pec, pic, ppc = conferenceProjectionStatistics(al, 4)
            print('array %d: b4 %s f4 %s' % (jj, b4, f4))
            print('  pec %s, pec5 %s, ppc %s, ppc5 %s' % (pec, pec5, ppc, ppc5))


#%%
if 0:
    rtable = createConferenceDesignsPageParetoTable(markup.page(), pareto_results, verbose=2, htmlsubdir=None)
    latextable = oapackage.array2latex(rtable)
    print(latextable)

#%%

if 0:
    from oapackage.oahelper import create_pareto_element
    from oaresearch.research_conference import createConferenceParetoElement, calculateConferencePareto

    presults = calculateConferencePareto(ll, N=None, k=None, verbose=1)


#%% Generate subpages for the designs


def conferenceSubPages(tag='conference', Nmax=40, Nstart=4, kmax=None,
                       verbose=1, specials={}, Nstep=2, NmaxPareto=40, cache=True, cache_tag='results_cachev6'):
    """ Generate a table with matrices

    Arguments:
        page (markup object)
        tag (str): type of data
        Nmax, Nstart, kmax (int)
        ta (str): text-align property
        verbose (int): verbosity level

    """

    Nrange = range(Nstart, Nmax + 1, Nstep)
    if kmax == -2:
        kmax = int(np.ceil(Nmax / 2) + 1)
    if kmax is None:
        kmax = np.max(Nrange) + 2
    krange = range(2, kmax)
    if verbose:
        print('conferenceSubPages: tag %s, Nmax %d, kmax %d' % (tag, Nmax, kmax))

    subpages = {}
    subpages[tag] = {}

    oapackage.mkdirc(os.path.join(outputdir, cache_tag))
    for ki, kk in enumerate(krange):
        for Ni, N in enumerate(Nrange):
            subpages[tag]['N%dk%d' % (N, kk)] = {}
            cachefile = os.path.join(outputdir, cache_tag, tag + '-' + 'N%dk%d' % (N, kk) + '.pickle')

            if cache and os.path.exists(cachefile):
                with open(cachefile, 'rb') as fid:
                    print('loading results from cachefile %s' % cachefile)
                    pareto_results, cfile = pickle.load(fid)
            else:

                # get arrays
                cfile, nn, mode = conferenceResultsFile(N, kk, outputdir, tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'], tagtype=['full', 'r', 'r'], verbose=1)

                if nn >= 10000 or N > NmaxPareto or mode != 'full':
                    continue

                ll = oapackage.readarrayfile(cfile)
                if verbose:
                    print('conferenceSubPages: generate %s N %d k %d: %d designs' % (tag, N, kk, nn))
                # calculate statistics
                presults, pareto = calculateConferencePareto(ll, N=N, k=kk, verbose=1)
                pareto_results = generateConferenceResults(presults, ll, ct=None, full=mode == 'full')
                pareto_results['arrayfile'] = cfile
                pareto_results['datadir'] = ''

                print('storing results in cachefile %s' % cachefile)

                with open(cachefile, 'wb') as fid:
                    pareto_results['presults'] = None
                    pickle.dump((pareto_results, cfile), fid)

            # create HTML page
            page = conferenceDesignsPage(pareto_results, verbose=1, makeheader=True, htmlsubdir=cdir)

            # write results
            htmlfile0 = os.path.basename(cfile).replace('.oa.gz', '.html').replace('.oa', '.html')
            htmlfile = os.path.join(cdir, htmlfile0)

            sx = subpages[tag]['N%dk%d' % (N, kk)]
            sx['htmlpage0'] = htmlfile0
            sx['htmlpage'] = htmlfile
            sx['pareto_results'] = pareto_results
            sx['arrayfile'] = cfile

            print('writing to %s' % htmlfile)
            with open(htmlfile, 'wt') as fid:
                fid.write(page)

    return subpages

#generated_subpages = conferenceSubPages(tag='cdesign', Nmax=40, Nstart=4, verbose=2, cache=False)
generated_subpages = conferenceSubPages(tag='cdesign', Nmax=40, Nstart=4, verbose=2, cache=True)

#%% Results table for latex

htmlsubdir = os.path.join(htmldir, 'conference')
for N in range(8,25,2):
    lst = oapackage.findfiles(htmlsubdir, 'conference-N%d.*pickle' % N)
    print('latex table: N %d: %d files' % (N, len(lst)))
    table = None

    kk = [oapackage.scanf.sscanf(file, 'conference-N%dk%d')[1] for file in lst]
    lst = [lst[idx] for idx in np.argsort(kk)]

    for file in (lst):
        r = pickle.load(open(os.path.join(htmlsubdir, file), 'rb'))

        ncolumns = r['ncolumns']
        rtable = r['rtable']
        if rtable.size==0:
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
    if len(lst)==0:
        print('no results for N=%d'  % N)
        continue

    offset_columns = [1, 2]
    for row in range(1, table.shape[0]):
        for col in offset_columns:
            table[row, col] = str(int(table[row, col]) + 1)
    latextable = oapackage.array2latex(table, hlines=[0], comment=['conference desgins N=%d' % (N), 'offset for indices is 1'])
    print(latextable)
    with open(os.path.join(htmlsubdir, 'conference-N%d-overview.tex' % (N, )), 'wt') as fid:
        fid.write(latextable)

#%%


def cdesignTag(N, kk, page, outputdir, tdstyle='', tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
               tagtype=['full', 'r', 'r'], verbose=1, ncache=None, subpage=None):
    """ Create html tag for oa page

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
        if 'full' not in ncache:
            ncache['full'] = {}
        ncache['full']['N%dk%d' % (N, kk)] = nn

    cfilex = oapackage.oahelper.checkOAfile(cfile)
    if cfilex is not None:
        cfilebase = os.path.basename(cfilex)
    else:
        cfilebase = None

    if page is not None:
        if subpage:
            hreflink = os.path.join('conference', subpage)
            print('hreflink: %s' % subpage)
        else:
            hreflink = 'conference/%s' % cfilebase

        txt, link = htmlTag(nn, kk, N, mode=mode,
                            href=hreflink, ncache=ncache, verbose=verbose >= 2)
        if verbose >= 2:
            print('cdesignTag: txt %s' % (txt,))
        if link:
            shutil.copyfile(cfilex, os.path.join(cdir, cfilebase))
        page.td(txt, style=tdstyle)
    else:
        if verbose >= 2:
            print(cfile)
    return cfile


#%% Testing
if 0:
    tag = 'dconferencej1'
    cdesignTag(N=N, kk=kk, page=page, outputdir=outputdir,
               tags=[tag, tag + '-r'], tagtype=['full', 'r'], verbose=2, ncache=ncache)

#%%


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
    s['cdesign'][24] = tt

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
#specialdata = {}

specialdataDC = specialdata

#%%


def DconferencePage(page, tag='dconference', Nmax=26, Nstart=4, kmax=None,
                    ta='left', verbose=1, specials={}, Nstep=2,
                    tableclass='conftable', tdstyle=None, subpages=None):
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
    subpage = None
    for ki, kk in enumerate(krange):
        page.tr()
        page.td('%s' % str(kk), style=tdstyle)
        for Ni, N in enumerate(Nrange):
            if subpages is not None:
                subpage = subpages.get('N%dk%d' % (N, kk), None)
                if subpage is not None and len(subpage) != 0:
                    subpage = subpage['htmlpage0']

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
                print('DconferencePage: tag %s: %d %d: subpage %s' % (tag, N, kk, subpage))

            if tag == 'cdesign':
                cdesignTag(N, kk, page, outputdir, tdstyle,
                           ['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'],
                           ['full', 'r', 'r'], subpage=subpage)
            else:
                cdesignTag(N, kk, page, outputdir, tdstyle,
                           [tag, tag + '-r'], ['full', 'r'], ncache=ncache)

        page.tr.close()
    page.table.close()
    page.p.close()

    return {'ncache': ncache}

#%% Get even-odd designs
if 0:
    Nrange = range(0, 82, 2)  # full range
    #Nrange=range(44, 45, 2)
    #Nrange=range(4, 48, 2)
    Nrange = range(74, 82, 2)
    tag = 'dconferencej1j3'
    for Ni, N in enumerate(Nrange):
        kmax = int(np.ceil(N / 2) + 1)
        krange = range(2, kmax)
        for ki, kk in enumerate(krange):

            cfile = cdesignTag(N, kk, page=None, outputdir=outputdir, tags=[
                               tag, tag + '-r'], tagtype=['full', 'r'])
            na = oapackage.nArrayFile(cfile)

            eolist = []
            if na > 100000:
                af = oapackage.arrayfile_t(cfile)
                for ii in range(na):
                    if ii % (200 * 1e3) == 0 or ii == na - 1:
                        print('N %d, k %d: %d/%d' % (N, kk, ii, af.narrays))
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

                # fixme: make streaming
                eolist = [
                    al for al in ll if not oapackage.isConferenceFoldover(al)]

            cfileout = cfile.replace(tag, tag + '-eo')
            print('  out: %s: %d -> %d' % (cfileout, na, len(eolist)))
            if 1:
                oapackage.writearrayfile(
                    cfileout, eolist, oapackage.ABINARY_DIFF, N, kk)
                xfile = cfileout + '.gz'
                if os.path.exists(xfile):
                    print('removing file %s' % (xfile))
                    os.remove(xfile)
                if 1:
                    if len(eolist) > 100:
                        cmd = 'gzip -f %s' % cfileout
                        os.system(cmd)

#DconferencePage(page, tag='dconferencej1j3-eo', Nstart=44, Nmax=45, kmax=28)


#%%


if generate_webpage:
    tdstyle = 'text-align:left; padding-left: .6em; margin-right: 0em; padding-right: .6em; margin-left: 0px;'

    citation = oaresearch.research.citation('cenumeration', style='brief')
    page = oapackage.markup.page()

    if not htemplate:
        page.init(title="Conference designs",
                  css=('oastyle.css'),
                  lang='en',
                  header="", htmlattrs=dict({'xmlns': 'http://www.w3.org/1999/xhtml'}),
                  # doctype=markup.doctype.strict,
                  doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"   "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">',
                  metainfo=({'keywords': 'conference desigms', 'robots': 'index, follow',
                             'description': 'Complete Enumeration of conference designs'}),
                  footer="<!-- End of page -->")

    page.h1('Non-isomorphic conference designs', style='margin-bottom: 0.1em;')
    page.p('<em style="color:darkblue;">pieter.eendebak@gmail.com</em>', style='margin-top: 0px;')
    ss = 'On this page we present numbers of isomorphism classes for conference designs with a specified number of runs (or rows) and factors (or columns). For all the cases a set of representatives for the isomorphism classes is available. The algorithm used to generate the results is described in %s. ' % citation
    #ss += 'For all the cases a set of representatives for the isomorphism classes in LMC-0 form is available.'
    page.p(ss)

    if 1:
        page.h2('Single conference designs')
        page.p()
        page.add(
            'Conference designs are matrixes of size 2m x k with values 0, +1 or -1, where m &ge; 2 and k &le; m. Each column contains exactly one zero and each row conains at most one zero. The columns are orthogonal to each other.')
        ss = 'Square single conference designs are <a href="https://en.wikipedia.org/w/index.php?title=Conference_matrix&oldid=702633334">conference matrices</a>.'
        page.add(ss)
        page.p.close()

        page.p('A definitive screening design can be constructed by appending a conference design with both its fold-over and a row of zeroes.')
        DconferencePage(page, tag='cdesign', Nmax=41,
                        kmax=32, specials=specialdata, tdstyle=tdstyle,
                        tableclass='conftable conftable1', subpages=generated_subpages['cdesign'])

    if 1:
        page.h2('Double conference designs (DCDs)')
        page.p()
        page.add('Double conference designs are matrices of size 4m x k with values 0, +1 or -1, where m &ge;2 and k &le; 2m. Each column contains exactly two zeros and each row contains at most one zero. The columns are orthogonal to each other.')
        page.p.close()

        def subheader(s):
            page.h3(s)
        subheader('DCDs with level balance and orthogonal interaction columns.')
        page.p('In double conference designs with level balance and orthogonal interaction columns, the elements of each column sums to zero. In addition, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero too. A definitive screening design can be constructed by appending a double conference design with a row of zeroes.')
        #page.p('All J1 and J3 values are zero.')
        DconferencePage(page, tag='dconferencej1j3', specials=specialdata,
                        Nstart=4, Nmax=81, kmax=28, Nstep=4, tableclass='conftable conftable2')
        #DconferencePage(page, tag='dconferencej1j3', Nstart=58, Nmax=81, kmax=22, Nstep=4)
        page.p()
        page.add('Only the even-odd designs')
        page.p.close()
        DconferencePage(page, tag='dconferencej1j3-eo',
                        Nstart=4, Nmax=81, kmax=24, Nstep=4, specials=specialdata, tableclass='conftable conftable2')
        #DconferencePage(page, tag='dconferencej1j3-eo', Nstart=58, Nmax=81, kmax=20)

        if 1:
            subheader('DCDs with level balance')  # J1 values are zero
            page.p('In double conference designs with level balance, the elements of each column sums to zero')
            r = DconferencePage(page, tag='dconferencej1', kmax=16, tableclass='conftable conftable3')

        if 1:
            subheader('DCDs with orthogonality of interaction columns')  # J3=0
            page.p('In double conference designs with orthogonal interaction columns, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero.')
            DconferencePage(page, tag='dconferencej3', kmax=16, Nmax=24, tableclass='conftable conftable4')

        if 1:
            subheader('Plain DCDs')  # no restrictions on J1 and J3
            page.p('In double conference designs with orthogonal interaction columns, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero.')
            DconferencePage(page, tag='dconference2', Nmax=20, kmax=16, tableclass='conftable conftable5')

        page.h2('Weighing matrices')

        if 1:
            subheader('Plain DCDs (full isomorphism class)')  # no restrictions on J1 and J3
            page.p()
            page.add(
                'For the isomorphism class level-permutatations of the rows are allowed.')
            page.add(
                'We also allow multiple zeros in a row.')
            page.add('The square double conference designs of this type form a complete non-isomorphic set of weighing matrices of type W(N, N-2).')
            page.p.close()
            # page.h2(
            #    'Double conference matrices with full isomorphism class (no restrictions on J1 and J3)')
            # page.p()
            # page.add(
            #    'Double conference matrices are matrixes of size 2m x k. Each columns contains exactly two zeros and an arbitrary number of values +1 or -1. The columns are orthogonal to each other.')
            page.add('The square double conference matrices are a complete non-isomorphic set of weighing matrices of type W(N, N-2).')
            page.p.close()
            DconferencePage(page, tag='weighing2', kmax=18, Nmax=25)

    ss = 'If you make use of these results, please cite the paper %s.' % citation
    page.p(ss)
    if 1:
        ss = 'The square conference designs are <a href="https://en.wikipedia.org/w/index.php?title=Conference_matrix&oldid=702633334">conference matrices</a>.'
        ss += 'The algorithm used to generate the results is described in %s.' % citation
        page.p(ss)

    if 1:
        page.h1('Obsolete results')
        subheader('DCDs with level balance and orthogonal interaction columns.')
        page.p('These designs have N &equiv; 2 mod 4. In double conference designs with level balance and orthogonal interaction columns, the elements of each column sums to zero. In addition, for any set of three columns, the vector formed by the element-wise product of the columns has elements that sum to zero too. A definitive screening design can be constructed by appending a double conference design with a row of zeroes.')
        DconferencePage(page, tag='dconferencej1j3', Nstart=6, Nmax=81, kmax=28, Nstep=4)

        if 0:
            subheader('Double conference matrices (no restrictions on J1 and J3)')
            page.p()
            page.add(
                'Double conference matrices are matrixes of size 2m x k. Each columns contains exactly two zeros and an arbitrary number of values +1 or -1. The columns are orthogonal to each other.')
            page.add(
                'The square double conference matrices are weighing matrices of type W(N, N-2).')
            page.p.close()
            DconferencePage(page, tag='dconference2')

    import webbrowser
    if htemplate:
        hfile = os.path.join(htmldir, 'templates', 'conference.html')
    else:
        hfile = os.path.join(htmldir, 'conference.html')

    with open(hfile, 'w') as fid:
        _ = fid.write(str(page))
    webbrowser.open_new_tab(hfile)
