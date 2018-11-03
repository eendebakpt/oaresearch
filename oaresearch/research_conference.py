# -*- coding: utf-8 -*-
""" Module to generate D-optimal designs

For more information see: https://doi.org/10.1080/00401706.2016.1142903

Pieter Eendebak <pieter.eendebak@gmail.com>

"""

from __future__ import print_function

import os
import numpy as np
import time
import itertools

try:
    import matplotlib
    import matplotlib.pyplot as plt
except:
    pass

import oapackage
import oapackage.markup as markup
import oapackage.oahelper as oahelper
from oapackage.markup import oneliner as e

import oaresearch
from oaresearch.research import citation
import oaresearch.filetools

#%%

def momentMatrix(k):
   """ Return the moment matrix of a conference design """
   pk=int(1+0.5*k*(k+1)+k)
   M=np.zeros((pk,pk));
   M[0,0]=1;
   M[0,int(0.5*k*(k+1)+1):]=1/3;
   M[int(0.5*k*(k+1)+1):,0]=1/3;
   M[1:(k+1),1:(k+1)]=np.eye(k)*1./3;
   M[(k+1):int(0.5*k*(k+1)+1),(k+1):int(0.5*k*(k+1)+1)]=np.eye(int(0.5*k*(k-1)))/9;
   M[int(0.5*k*(k+1)+1):,int(0.5*k*(k+1)+1):]=np.eye(k)/5+(np.ones( (k,k) )-np.eye(k))/9;
   return M

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
            #print('%d %d -> %d'  %(ii, jj,idx) )
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

    # print(al.Jcharacteristics(4))
    #X2 = al.getModelMatrix(2)[:, (1 + al.n_columns):]
    X2 = conferenceSecondOrder(al, False)
    rank = np.linalg.matrix_rank(X2)
    # oapackage.makearraylink(X2).showarray()

    X2 = conferenceSecondOrder(al, True)
    rankq = np.linalg.matrix_rank(X2)
    # oapackage.makearraylink(X2).showarray()

    if verbose:
        x = np.array(al)
        x2 = x.copy()
        x2 *= -1
        # x2[:,0]*=-1
        foldover = np.vstack((x, x2))
        foldover = oapackage.makearraylink(foldover)
        #X2 = al.getModelMatrix(2)[:, (1 + al.n_columns):]
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


def test_makePareto():
    pass


class designResults:
    pass


def calculateConferencePareto(ll, N=None, k=None, verbose=1):
    if verbose:
        print('analysing %d arrays ' % (len(ll), ))
        print('  criteria: rank, F4, B4')
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

    pareto = makePareto(presults, addFoldover=True)

    presults.pareto_indices = pareto.allindices()
    presults.nclasses = pareto.number()
    presults.npareto = pareto.numberindices()

    # print(al.GWLP()[0:5])
    return presults, pareto


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
    rr = {'type': 'conference designs', 'arrayfile': None, 'presults': presults}
    rr['full'] = full
    rr['N'] = presults.N
    rr['k'] = presults.k
    rr['idstr'] = 'cdesign-%d-%d' % (rr['N'], rr['k'])
    if ct is not None:
        rr['ctidstr'] = ct.idstr()
        assert(ct.N == presults.N)
        assert(ct.ncols == presults.k)

    rr['narrays'] = len(ll)
    if len(ll) == 0:
        pass
    else:
        rr['minrankX2'] = np.min(presults.ranks)
        rr['maxrankX2'] = np.max(presults.ranks)
        rr['minrankX2q'] = np.min(presults.ranksX2Q)
        rr['maxrankX2q'] = np.max(presults.ranksX2Q)

        rr['maxB4'] = np.max(presults.b4s)
        rr['minB4'] = np.min(presults.b4s)

        rr['minF4'] = presults.f4s[oapackage.oahelper.sortrows(np.abs(presults.f4s))[0]]
        rr['maxF4'] = presults.f4s[oapackage.oahelper.sortrows(np.abs(presults.f4s))[-1]]
    if hasattr(presults, 'pareto_indices'):
        rr['pareto_designs'] = [ll[ii] for ii in presults.pareto_indices]

    return rr

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


def conferenceDesignsPage(rr, verbose=1, makeheader=True, htmlsubdir=None):
    """ Generate html page for class conference desgins """
    #idstr = rr['ctidstr']
    #idstr = '%s' % (ct.idstr(),)
    N = rr['N']
    k = rr['k']

    if makeheader:
        pext = 'html'
    else:
        pext = 'html.template'

    xstr = 'C(%d, %d)' % (N, k)
    xstrplain = xstr

    if verbose:
        print('conferenceDesignsPage: start of generation')

    presults = rr.get('presults', None)

    page = markup.page()
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

    full_results = rr.get('full')
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

    narrays = rr['narrays']
    simpleRow('Number of non-isomorphic designs', str(rr['narrays']))

    if narrays > 0:
        simpleRow('Minimum/Maximum rank of X2', '%d/%d' % (rr['minrankX2'], rr['maxrankX2']))
        simpleRow('Minimum/Maximum rank of X2 with quadratics', '%d/%d' % (rr['minrankX2q'], rr['maxrankX2q']))
        simpleRow('Minimum B4', '%.4f' % rr['minB4'])
        simpleRow('Maximum B4', '%.4f' % rr['maxB4'])
        simpleRow('Minimum F4', str(rr['minF4']))
        simpleRow('Maximum F4', str(rr['maxF4']))

    havearrayfile = 0
    if 'arrayfile' in list(rr.keys()):
        if not rr['arrayfile'] is None:
            havearrayfile = 1

    if havearrayfile:
        iafile = rr['arrayfile']
        outfile0 = rr['idstr'] + '.oa'

        na = oapackage.nArrayFile(iafile)
        if verbose >= 2:
            print('conferenceDesignsPage: read arrayfile %s: na %d' % (iafile, na))

        if na < 5000 and na >= 0 and 1:
            if htmlsubdir is None:
                raise Exception('need html subdirectory to copy .oa file')
            outfilefinal = oaresearch.filetools.copyOAfile(iafile, htmlsubdir, outfile0, convert='T', zipfile=None, verbose=1, cache=0)

            if rr.get('full', False):
                htag = oaresearch.htmltools.formatArrayHyperlink(
                    'all arrays', outfile0, iafile)

                rr['datafilestr'] = 'all arrays'
            else:
                na = oapackage.nArrayFile(os.path.join(htmlsubdir, outfilefinal))
                htag = oaresearch.htmltools.formatArrayHyperlink(
                    'all arrays', outfile0, iafile)
                rr['datafilestr'] = '%d array(s)' % na

            simpleRow('Data', htag)

        else:
            if verbose:
                print('conferenceDesignsPage: no datafile (na %d)' % na)
            rr['datafilestrfile'] = None
            rr['datafilestr'] = '-'

    if 'totaltime' in list(rr.keys()):
        simpleRow('Processing time', str(rr['totaltime']) + 's')
    page.table.close()

    if narrays > 0 and full_results:
        add_extra = True
        print('do statistics2htmltable')

        header = ['Array index', 'Rank X2', 'F4', 'B4']
        if add_extra:
            header += ['Rank X2 with q']

        rtable = np.zeros((1 + len(presults.pareto_indices), len(header)), dtype='|U208')
        rtable[:] = ' '
        for ii, h in enumerate(header):
            rtable[0, ii] = header[ii]

        pidx = presults.pareto_indices  # pareto.allindices()
        for ii, idx in enumerate(pidx):
            rtable[ii + 1, 0:4] = ['%d' % idx, str(presults.ranks[idx]), str(presults.f4s[idx]), '%.4f' % (presults.b4s[idx])]
            if add_extra:
                rtable[ii + 1, -1:] = [str(presults.ranksX2Q[idx])]

        subpage = oaresearch.research.array2html(rtable, header=1, tablestyle='border-collapse: collapse;',
                                                 trclass='', tdstyle='padding-right:1em;', trstyle='', thstyle='text-align:left; padding-right: 1em;', comment=None)
        page.br(clear='both')
        page.h2('Pareto optimal designs')
        page.p()
        page.add('There are %d Pareto optimal designs in %d classes.' % (presults.npareto, presults.nclasses))
        page.add('Pareto optimality is according to rank, F4 and B4 (the other statistics are ignored).')
        page.p.close()
        if rr.get('pareto_designs', None) is not None:
            pdesigns = rr.get('pareto_designs', None)

        pfile0 = outfile0 = rr['idstr'] + '-pareto.oa'

        if htmlsubdir is not None:
            pfile = os.path.join(htmlsubdir, pfile0)
            oapackage.writearrayfile(pfile, pdesigns)
            page.p('All %s' % e.a('Pareto optimal designs', href=pfile0) + '.')

        page.add(str(subpage))
        page.br(clear='both')


    localtime = time.asctime(time.localtime(time.time()))
    dstr = str(localtime)

    page.p('<br/>\n')  
    page.p('Page generated on %s.' % dstr)

    pstr = str(page).replace('<meta content="charset=utf-8" name="text/html" />',
                             '<meta http-equiv="Content-Type" content="charset=utf-8" name="text/html" />')

    return pstr
