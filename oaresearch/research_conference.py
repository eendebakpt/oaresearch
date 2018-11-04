# -*- coding: utf-8 -*-
""" Module to generate D-optimal designs

For more information see: https://doi.org/10.1080/00401706.2016.1142903

Pieter Eendebak <pieter.eendebak@gmail.com>

"""

from __future__ import print_function

import os
import sys
import numpy as np
import time
import itertools

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

#%%

from oapackage.conference import momentMatrix, modelStatistics, conferenceProjectionStatistics


def generateConference(N, kmax=None, verbose=1, diagc=False, nmax=None, selectmethod='random', tag='cdesign', outputdir=None):
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
    if not nmax is None:
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
            print('generateConference: extcol %d: %d designs' % (extcol, len(LL[extcol - 1])) )
            sys.stdout.flush()
        LL[extcol] = oapackage.extend_conference(
            LL[extcol - 1], ctype, verbose=verbose >= 2)

        LL[extcol] = oapackage.selectConferenceIsomorpismClasses(LL[extcol], 1)

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
                    # mixed
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

   
#%%
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
    b4 = np.sum(np.array(j4)**2) / N**2

    ncols=al.n_columns
    X2 = oapackage.array2modelmatrix(al, 'i')[:,(1+ncols):]
    rank = np.linalg.matrix_rank(X2)

    X2 = oapackage.array2modelmatrix(al, 'i')
    rankq = np.linalg.matrix_rank(X2)

    if verbose:
        x = np.array(al)
        x2 = x.copy()
        x2 *= -1
        foldover = np.vstack((x, x2))
        foldover = oapackage.makearraylink(foldover)
        X2 = conferenceSecondOrder(foldover, False)
        rankfoldover = np.linalg.matrix_rank(X2)

    if verbose:
        print('f4: %s' % (f4,))
        print('j4: %s' % (j4,))
        print('rank X2: %s' % (rank,))
        print('rank X2+quadratics: %s' % (rankq,))
        print('rank of foldover X2: %s' % (rankfoldover,))
    return [f4, b4, rank, rankq]


def test_confJ4():
    al = oapackage.exampleArray(18)
    J = conferenceJ4(al)
    assert(np.sum(np.abs(np.array(J)) == 12) == 1)
    assert(np.sum(np.abs(np.array(J)) == 0) == 23)

from oapackage.oahelper import create_pareto_element

def createConferenceParetoElement(al, addFoldover=True, addProjectionStatistics=True, pareto=None):
    """ Create Pareto element from conference design """
    rr = oaresearch.research_conference.conferenceStatistics(al, verbose=0)
    [f4, b4, rankinteraction, ranksecondorder] = rr[0:4]
    f4minus = [-float(x) for x in f4]
    values = [ list(f4minus), [-float(b4)], [float(ranksecondorder)] ]
    data ={'F4': f4, 'B4': b4, 'ranksecondorder': ranksecondorder}
    
    for kk in [4,5]:
        pec, pic, ppc = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=kk)

        if addProjectionStatistics:
            values += [ [pec] ]
            values += [ [pic] ]
            values += [ [ppc] ]
        data['pec%d' % kk] = pec
        data['pic%d' % kk] = pic
        data['ppc%d' % kk] = ppc


    if addFoldover:
            foldover=oapackage.isConferenceFoldover(al)
            values += [[int(foldover)], [int(not foldover)]]

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




class designResults:
    pass

    def add_value(self, tag, value):
        mintag=tag+'_min'
        maxtag=tag+'_max'
        if not hasattr(self, mintag):
                setattr(self, mintag, value)
        if not hasattr(self, maxtag):
                setattr(self, maxtag, value)
        setattr(self, mintag, min(value, getattr(self, mintag)))
        setattr(self, maxtag, max(value, getattr(self, maxtag)))

def calculateConferencePareto(ll, N=None, k=None, verbose=1, add_data = True):
    """ Calculate Pareto optimal designs from a list of designs """
    if verbose:
        print('calculateConferencePareto: analysing %d arrays ' % (len(ll), ))
    pareto = oapackage.ParetoMultiDoubleLong()
    
    presults = designResults()
    for ii, al in  enumerate(ll):
        oapackage.oahelper.tprint('calculateConferencePareto: N %s column %s: array %d/%d: %s' % (str(N), str(k), ii, len(ll), str(pareto).strip() ) , dt=2)
        pareto_element, data = createConferenceParetoElement(al, addFoldover=True, pareto=pareto)
        
        pareto.addvalue(pareto_element, ii)

        if add_data:
            for tag in ['ranksecondorder', 'B4', 'pec4', 'pic4']:
                presults.add_value(tag, data[tag])
    presults.N = N
    presults.ncolumns = k
    if len(ll) > 0:
        presults.N = ll[0].n_rows
        presults.ncolumns = ll[0].n_columns
        
    pareto.show()
    presults.pareto_indices = pareto.allindices()
    presults.nclasses = pareto.number()
    presults.npareto = pareto.numberindices()

    presults.pareto_designs = [ll[ii] for ii in presults.pareto_indices]
    presults.pareto_data=[]
    for ii, al in  enumerate(presults.pareto_designs):
        pareto_element, data = createConferenceParetoElement(al, addFoldover=True)
        presults.pareto_data.append(data)

    return presults, pareto

def test_calculateConferencePareto():
    ll=[oapackage.exampleArray(idx) for idx in [45,46,47,45]]
    presults, _ = calculateConferencePareto(ll, N=None, k=None, verbose=1, add_data = True)

if __name__=='__main__':
    test_calculateConferencePareto()
    
def calculateConferenceParetoX(ll, N=None, k=None, verbose=1):
    if verbose:
        print('analysing %d arrays ' % (len(ll), ))
        print('  criteria: rank, F4, B4, ...')
        print('\n')

    presults = designResults()
    presults.N = N
    presults.k = k
    if len(ll) > 0:
        presults.N = ll[0].n_rows
        presults.k = ll[0].n_columns

    presults.ranks = []
    presults.ranksX2Q = []
    presults.b4s = []
    presults.f4s = []
    presults.foldover = []
    presults.pec4=[]
    presults.pic4=[]
    presults.ppc4=[]
    for ii, al in enumerate(ll):
        if verbose:
            oapackage.tprint('processing array %d/%d' % (ii, len(ll)), dt=5)

        rr = conferenceStatistics(al, verbose=0)
        [f4, b4, rank, rankq] = rr[0:4]
        N = al.n_rows

        presults.b4s += [b4]
        presults.f4s += [f4]
        presults.ranks += [rank]
        presults.ranksX2Q += [rankq]
        presults.foldover += [oapackage.isConferenceFoldover(al)]
        pec4, pic4, ppc4 = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=4)
        pec5, pic5, ppc5 = oapackage.conference.conferenceProjectionStatistics(al, ncolumns=4)

        presults.pec4+=[pec4]
        presults.pic4+=[pic4]
        presults.ppc4+=[ppc4]
        
    pareto = makePareto(presults, addFoldover=True)

    presults.pareto_indices = pareto.allindices()
    presults.nclasses = pareto.number()
    presults.npareto = pareto.numberindices()

    return presults, pareto

def showMaxZ(LL):
    """ For a list of generated designs show the maximum zero position """
    N = LL[3][0].n_rows

    for ii, L in enumerate(LL):
        k = ii + 1
        s = [oapackage.maxz(al) for al in L]
        mm, _ = np.histogram(s, range(N + 1))
        print('%d cols: maxz seq %s' % (k, list(mm)) )
        
def conferenceResultsFile(N, kk, outputdir, tags=['cdesign', 'cdesign-diagonal', 'cdesign-diagonal-r'], tagtype=['full', 'r', 'r'], verbose=1):
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
    pareto_results = {'type': 'conference designs', 'arrayfile': None, 'presults': presults}
    pareto_results['full'] = full
    pareto_results['full_results'] = full
    pareto_results['N'] = presults.N
    pareto_results['ncolumns'] = presults.ncolumns
    pareto_results['idstr'] = 'cdesign-%d-%d' % (pareto_results['N'], pareto_results['ncolumns'])
    if ct is not None:
        pareto_results['ctidstr'] = ct.idstr()
        assert(ct.N == presults.N)
        assert(ct.ncols == presults.k)

    pareto_results['narrays'] = len(ll)
    if len(ll) == 0:
        pass
    else:
        pareto_results['minrankX2q'] =(presults.ranksecondorder_min)
        pareto_results['maxrankX2q'] = (presults.ranksecondorder_max)
        #rr['minrankX2q'] = np.min(presults.ranksX2Q)
        #rr['maxrankX2q'] = np.max(presults.ranksX2Q)

        pareto_results['maxB4'] = (presults.B4_min)
        pareto_results['minB4'] = (presults.B4_max)

        #rr['minF4'] = presults.f4s[oapackage.oahelper.sortrows(np.abs(presults.f4s))[0]]
        #rr['maxF4'] = presults.f4s[oapackage.oahelper.sortrows(np.abs(presults.f4s))[-1]]
    #if hasattr(presults, 'pareto_indices'):
    pareto_results['pareto_designs'] = getattr(presults, 'pareto_designs', None) # = [ll[ii] for ii in presults.pareto_indices]
    pareto_results['pareto_designs'] = [np.array(array) for array in pareto_results['pareto_designs']]
    pareto_results['pareto_indices'] = presults.pareto_indices
    for tag in ['npareto', 'nclasses', 'pareto_data']:
        pareto_results[tag]=getattr(presults, tag)
    return pareto_results

#%% Webpage generation


def nprevzero(N, k, ncache):
    """ Return true if any previous result was zero """
    for ix in range(k - 1, 2, -1):
        # print(ix)
        p = ncache['full'].get('N%dk%d' % (N, ix), -1)
        #p = ncache['full'].get('N%dk%d' % (N, ix), -1)
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

    """
    link = False
    if nn >= 0:
        if mode == 'full':
            txt = '%d' % nn
        else:
            if nn == 0:
                txt = '?'
            else:
                txt = '&ge; %d' % nn
        if nn < 10000 and nn > 0:
            if href is not None:
                ss = e.a(txt, href=href, style='text-decoration: none;')
                link = True
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
    return txt, link


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


def createConferenceDesignsPageHeader(page, makeheader, conference_class, ncolumns):
    
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

    page.h1('Conference designs %s ' % xstr)
    oap = e.a('Orthogonal Array package', href='../../software.html')
    pstr = 'This page contains information about conference designs. '
    pstr += 'The results have been generated with the %s package.' % oap
    pstr += ' If you use these data, please cite the paper ' + \
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
        #simpleRow('Minimum/Maximum rank of X2', '%d/%d' % (rr['minrankX2'], rr['maxrankX2']))
        simpleRow('Minimum/Maximum rank of X2 with quadratics', '%d/%d' % (pareto_results['minrankX2q'], pareto_results['maxrankX2q']))
        simpleRow('Minimum B4', '%.4f' % pareto_results['minB4'])
        simpleRow('Maximum B4', '%.4f' % pareto_results['maxB4'])
        #simpleRow('Minimum F4', str(rr['minF4']))
        #simpleRow('Maximum F4', str(rr['maxF4']))
    if 'totaltime' in list(pareto_results.keys()):
        simpleRow('Processing time', str(pareto_results['totaltime']) + 's')
        
    simpleRow('Data', pareto_results['datafile_tag'])

    page.table.close()
    
    
def createConferenceDesignsPageLoadDesignsFile( pareto_results, htmlsubdir=None, verbose=1):

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
            outfilefinal = oaresearch.filetools.copyOAfile(iafile, htmlsubdir, outfile0, convert='T', zipfile=None, verbose=1, cache=0)

            if pareto_results.get('full', False):
                htag = oaresearch.htmltools.formatArrayHyperlink(
                    'all arrays', outfile0, iafile)

                pareto_results['datafilestr'] = 'all arrays'
                pareto_results['datafile_tag'] = htag
            else:
                na = oapackage.nArrayFile(os.path.join(htmlsubdir, outfilefinal))
                htag = oaresearch.htmltools.formatArrayHyperlink(
                    'all arrays', outfile0, iafile)
                pareto_results['datafilestr'] = '%d array(s)' % na
                pareto_results['datafile_tag'] = htag


        else:
            if verbose:
                print('conferenceDesignsPage: no datafile (na %d)' % na)
            pareto_results['datafilestrfile'] = None
            pareto_results['datafilestr'] = '-'

def createConferenceDesignsPageParetoTable(page, pareto_results, verbose=0, htmlsubdir = None):

    if verbose:
        print('createConferenceDesignsPageParetoTable: start')
        
    pareto_indices = pareto_results['pareto_indices']
    pareto_data = pareto_results['pareto_data']
        
    if pareto_results['narrays'] > 0 and pareto_results.get('full_results'):
        add_extra = True
        print('do statistics2htmltable')

        header = ['Array index', 'Rank X2', 'F4', 'B4']
        if add_extra:
                for kk in [4,5]:
                    for tag  in ['PEC', 'PIC', 'PPC']:
                        header+=[tag +'%d' % kk]
                        
        rtable = np.zeros((1 + len(pareto_results['pareto_indices']), len(header)), dtype='|U208')
        rtable[:] = ' '
        for ii, h in enumerate(header):
            rtable[0, ii] = header[ii]

        for pareto_idx, list_idx in enumerate(pareto_indices):
            rtable[pareto_idx + 1, 0:4] = ['%d' % pareto_idx,  str(pareto_data[pareto_idx]['ranksecondorder']), str(pareto_data[pareto_idx]['F4']), '%.4f' % ((pareto_data[pareto_idx]['B4']))]
            if add_extra:
                column_offset=4
                for kk in [4,5]:
                    for tag  in ['pec', 'pic', 'ppc']:
                        rtable[pareto_idx + 1, column_offset] = '%.4f' % (pareto_data[pareto_idx][tag +'%d' % kk])
                        column_offset=column_offset+1
                        
        subpage = oaresearch.research.array2html(rtable, header=1, tablestyle='border-collapse: collapse;',
                                                 trclass='', tdstyle='padding-right:1em;', trstyle='', thstyle='text-align:left; padding-right: 1em;', comment=None)
        page.br(clear='both')
        page.h2('Pareto optimal designs')
        page.p()
        if pareto_results['npareto']==1:
            page.add('There is %d Pareto optimal designs in %d classes.' % (pareto_results['npareto'], pareto_results['nclasses']))
        else:
            page.add('There are %d Pareto optimal designs in %d classes.' % (pareto_results['npareto'], pareto_results['nclasses']))
        page.add('Pareto optimality is according to rank, F4 and B4, PECk, PICk and PPCk (the other statistics are ignored).')
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

            
def conferenceDesignsPage(pareto_results, verbose=1, makeheader=True, htmlsubdir=None):
    """ Generate html page for class conference designs """

    N = pareto_results['N']
    ncolumns = pareto_results['ncolumns']
   

    if makeheader:
        pext = 'html'
    else:
        pext = 'html.template'

    conference_class=oapackage.conference_t(N, ncolumns, 0)


    if verbose:
        print('conferenceDesignsPage: start of generation')

    page = markup.page()
    createConferenceDesignsPageHeader(page, makeheader=makeheader, conference_class=conference_class, ncolumns=ncolumns)
    createConferenceDesignsPageLoadDesignsFile( pareto_results, htmlsubdir=htmlsubdir)
    createConferenceDesignsPageResultsTable(page, pareto_results, verbose=verbose)
    createConferenceDesignsPageParetoTable(page, pareto_results, verbose=verbose, htmlsubdir=htmlsubdir)

    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)

    page.p('<br/>\n')  
    page.p('Page generated on %s.' % dstr)

    pstr = str(page).replace('<meta content="charset=utf-8" name="text/html" />',
                             '<meta http-equiv="Content-Type" content="charset=utf-8" name="text/html" />')

    return pstr
