# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

#%%
import numpy as np
import unittest
import collections

import oapackage
import oaresearch.research_conference

#%%


class TestResearchConference(unittest.TestCase):

    def test_createConferenceParetoElement(self):
        al1 = oapackage.exampleArray(49)
        pareto1, data1 = oaresearch.research_conference.createConferenceParetoElement(al1)
        al2 = oapackage.exampleArray(50)
        pareto2, data2 = oaresearch.research_conference.createConferenceParetoElement(al2)

        self.assertEqual(data1, collections.OrderedDict([('ranksecondorder', 20),
                                                         ('rankinteraction', 19),
                                                         ('F4', (0, 0, 4, 46, 20)),
                                                         ('B4', 2.48),
                                                         ('PEC4', 1.0),
                                                         ('PEC5', 1.0),
                                                         ('PIC4', 17.105),
                                                         ('PIC5', 16.537),
                                                         ('PPC4', 1.789),
                                                         ('PPC5', 1.213),
                                                         ('foldover', 0),
                                                         ('notfoldover', 1)]))

        # numerial instabilities can cause rounding errors in the PPC5. this should be solved by rounding in test_createConferenceParetoElement
        self.assertEqual(data1['PPC5'], data2['PPC5'])

    def test_generateConference(self):
        LL = oaresearch.research_conference.generateConference(8, verbose=0)
        self.assertEqual([len(l) for l in LL], [0, 1, 1, 2, 1, 1, 1, 1])

    def test_calculateConferencePareto(self):
        arrays = [oapackage.exampleArray(idx, 0) for idx in [45, 46, 47, 48]]

        presults, pareto = oaresearch.research_conference.calculateConferencePareto(arrays, N=None, k=None, verbose=1)
        self.assertEqual(presults['nclasses'], 2)
        self.assertEqual(presults['pareto_indices'], (0, 3))


    def test_conferenceStatistics(self):
        array = oapackage.exampleArray(51)
        expected = [(0, 1, 0), 0.1111111111111111, 6, 9]
        results = oaresearch.research_conference.conferenceStatistics(array)
        self.assertEqual(expected, results)

        array = oapackage.exampleArray(52)
        expected = [(0, 0, 1), 0.0, 6, 10]
        results = oaresearch.research_conference.conferenceStatistics(array)
        self.assertEqual(expected, results)

    def test_conference_statistics(self):
        arrays = [oapackage.exampleArray(idx, 0) for idx in [45, 46, 47, 48]]

        j4s = [oaresearch.research_conference.conferenceJ4(al) for al in arrays]
        self.assertEqual([np.sum(j4) for j4 in j4s], [52, 36, 36, 144])

        N = arrays[0].n_rows
        b4s = [np.sum(np.array(j4)**2) / N**2 for j4 in j4s]
        self.assertEqual(b4s, [3.16, 3.72, 3.4, 6.0])

        F4_values = [al.FvaluesConference(number_of_columns=4) for al in arrays]
        self.assertEqual(F4_values, [(0, 1, 4, 54, 11), (0, 1, 8, 52, 9), (0, 2, 4, 51, 13), (0, 0, 27, 42, 1)])


if __name__ == '__main__':
    unittest.main()
