# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

#%%
import numpy as np
import unittest

import oapackage
import oaresearch.research_conference

#%%


class TestResearchConference(unittest.TestCase):

    def test_generateConference(self):
        LL=oaresearch.research_conference.generateConference(8)
        self.assertEqual([len(l) for l in LL], [0, 1, 1, 2, 1, 1, 1, 1])
        
    def test_calculateConferencePareto(self):
        arrays=[oapackage.exampleArray(idx,1 ) for idx in [45,46,47,48]]

        presults, pareto = oaresearch.research_conference.calculateConferencePareto(arrays, N=None, k=None, verbose=1)
        self.assertEqual(presults['nclasses'], 2)
        self.assertEqual(presults['pareto_indices'], (0,3) )
        
    def test_conference_statistics(self):
        arrays=[oapackage.exampleArray(idx,1 ) for idx in [45,46,47,48]]
        
        j4s=[oaresearch.research_conference.conferenceJ4(al) for al in arrays]
        self.assertEqual([np.sum(j4) for j4 in j4s], [52, 36, 36, 144])
        
        N=arrays[0].n_rows
        b4s = [np.sum(np.array(j4)**2) / N**2 for j4 in j4s]
        self.assertEqual(b4s, [3.16, 3.72, 3.4, 6.0])
        
        F4_values = [al.FvaluesConference(jj=4) for al in arrays]
        self.assertEqual(F4_values, [(0, 1, 4, 54, 11), (0, 1, 8, 52, 9), (0, 2, 4, 51, 13), (0, 0, 27, 42, 1)])        


if __name__ == '__main__':
    unittest.main()
