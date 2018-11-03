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
        
        N=arrays[0].n_rows
        b4s = [np.sum(np.array(j4)**2) / N**2 for j4 in j4s]
        self.assertEqual(b4s, [3.16, 3.72, 3.4, 6.0])
        
        F4_values = [al.FvaluesConference(jj=4) for al in arrays]
        self.assertEqual(F4_values, [(0, 1, 4, 54, 11), (0, 1, 8, 52, 9), (0, 2, 4, 51, 13), (0, 0, 27, 42, 1)])        

    def test_conference_projection_statistics(self):
        arrays=[oapackage.exampleArray(idx,1 ) for idx in [45,46,47,48]]
        statistics = [ oaresearch.research_conference.conferenceProjectionStatistics(array) for array in arrays ]
        self.assertEqual(statistics, [(0.9857142857142858, 16.831420420862223, 1.7523727421346018), (0.9857142857142858, 16.782711360159983, 1.7515927378529395), (0.9714285714285714, 16.594006884359878, 1.7298655460036123), (1.0, 16.80031370680457, 1.7750342960753236)])
    
    def test_leftDivide(self):
        A = np.array([[1,2],[3,4]])
        B = np.array([[51,-2],[3,4]])
        C= oaresearch.research_conference.leftDivide(A, B)
        np.testing.assert_array_almost_equal(C,np.array([[-99.,   8.], [ 75.,  -5.]]))

if __name__ == '__main__':
    unittest.main()
