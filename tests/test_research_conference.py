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

    def test_momentMatrix(self):
        M1=np.array([[1.        , 0.        , 1./3],
       [0.        , 1./3, 0.        ],
       [1./3, 0.        , 0.2       ]])
        M=oaresearch.research_conference.momentMatrix(1)
        np.testing.assert_array_equal(M,M1)
        M=oaresearch.research_conference.momentMatrix(4)
        self.assertEqual(M.shape, (15,15))
        
    def test_conference_statistics(self):
        arrays=[oapackage.exampleArray(idx,1 ) for idx in [45,46,47,48]]
        
        j4s=[oaresearch.research_conference.conferenceJ4(al) for al in arrays]
        self.assertEqual([np.sum(j4) for j4 in j4s], [52, 36, 36, 144])
        F4_values = [al.FvaluesConference(jj=4) for al in arrays]
        self.assertEqual(F4_values, [(0, 1, 4, 54, 11), (0, 1, 8, 52, 9), (0, 2, 4, 51, 13), (0, 0, 27, 42, 1)])        

if __name__ == '__main__':
    unittest.main()
