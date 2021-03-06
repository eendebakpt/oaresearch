# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

# %%
import numpy as np
import unittest
import collections
import unittest.mock as mock
from unittest.mock import patch
from collections import OrderedDict
import io

import oapackage
import oaresearch.research_conference
from oaresearch.research_conference import (conference_design_has_extensions, conference_design_has_maximal_extension,
                                            maximal_extension_size, createConferenceParetoElement, DesignStack)

# %%


class TestFullExtensions(unittest.TestCase):

    def test_conference_design_has_maximal_extension(self):
        conference_class = oapackage.conference_t(12, 12, 0)
        conference_designs = [conference_class.create_root_three_columns()]
        all_data = {3: conference_designs}

        for ii, ncols in enumerate(range(4, conference_class.ncols + 1)):
            arrays = oapackage.extend_conference(all_data[ncols - 1], conference_class, verbose=0)
            arrays = oapackage.selectConferenceIsomorpismClasses(arrays, verbose=0)
            all_data[ncols] = arrays

        all_data_nauty = {}
        sort_idx={}
        data_nauty_sorted={}
        for k in all_data:
            all_data_nauty[k] = [oapackage.reduceConference(d) for d in all_data[k]]
            sort_idx[k] = np.argsort(np.array(all_data_nauty[k], dtype=object))
            data_nauty_sorted[k]= [ all_data_nauty[k][idx] for idx in sort_idx[k]]
        design_stack = DesignStack(all_data, all_data_nauty, sort_idx, data_nauty_sorted)

        number_of_classes = [len(value) for key, value in design_stack[0].items()]

        self.assertEqual(number_of_classes, [1, 3, 2, 5, 2, 2, 1, 1, 1, 1])

        conference_design_has_maximal_extension.design_stack = design_stack
        results = [conference_design_has_maximal_extension(array) for array in design_stack[0][4]]
        self.assertEqual(results, [True, True, False])

        results = [conference_design_has_maximal_extension(array) for array in design_stack[0][6]]
        self.assertEqual(results, [True, False, True, True, True])


    def test_maximal_extension_size(self):
        array = oapackage.exampleArray(40)
        m, extensions = maximal_extension_size(array, verbose=0)
        self.assertEqual(m, 14)
        self.assertEqual([array.md5() for aray in extensions], [
                         '045767c7cf6a005c06a6101dcf546ae0'])

        array = oapackage.exampleArray(50)
        m, extensions = maximal_extension_size(array, verbose=0)
        self.assertEqual(m, 12)
        self.assertEqual([array.md5() for aray in extensions], [
                         '72f829905c4fea1bf17a8d37686ac813'])


class TestCreateConferenceParetoElement(unittest.TestCase):

    def test_createConferenceParetoElement(self):
        array = oapackage.exampleArray(49)

        pareto_element, data = createConferenceParetoElement(array)
        pareto_data = [list(e.values) for e in list(pareto_element)]

        self.assertIsInstance(data, OrderedDict)
        self.assertEqual(pareto_data, [[20.0],
                                       [19.0],
                                       [-0.0, -0.0, -4.0, -46.0, -20.0],
                                       [-2.48],
                                       [1.0],
                                       [1.0],
                                       [17.105],
                                       [16.537],
                                       [1.789],
                                       [1.213],
                                       [0.0],
                                       [1.0]])

        self.assertDictEqual(data, OrderedDict([('ranksecondorder', 20),
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
        pareto_element, data = createConferenceParetoElement(
            array, addMaximumExtensionColumns=True)
        self.assertDictEqual(data, OrderedDict([('ranksecondorder', 20),
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
                                                ('notfoldover', 1),
                                                ('maximum_extension_size', 8)]))


class TestResearchConference(unittest.TestCase):

    def test_createConferenceParetoElement(self):
        al1 = oapackage.exampleArray(49)
        pareto1, data1 = oaresearch.research_conference.createConferenceParetoElement(
            al1)
        pareto1_data = [list(e.values) for e in list(pareto1)]
        self.assertEqual(pareto1_data[0], [20.])

        al2 = oapackage.exampleArray(50)
        pareto2, data2 = oaresearch.research_conference.createConferenceParetoElement(
            al2)
        pareto2_data = [list(e.values) for e in list(pareto1)]

        self.assertEqual(pareto2_data[0], [20.])

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

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            presults, pareto = oaresearch.research_conference.calculateConferencePareto(
                arrays, N=None, k=None, verbose=1)
            stdout = mock_stdout.getvalue()
            self.assertTrue(stdout.startswith('calculateConferencePareto: analysing 4 arrays'))
        self.assertEqual(presults['nclasses'], 2)
        self.assertEqual(presults['pareto_indices'], (0, 3))
        self.assertEqual(presults['pareto_data'][0], OrderedDict([('ranksecondorder', 20),
                                                                  ('rankinteraction', 19),
                                                                  ('F4', (0, 1,
                                                                          4, 54, 11)),
                                                                  ('B4', 3.16),
                                                                  ('PEC4', 0.986),
                                                                  ('PEC5', 0.929),
                                                                  ('PIC4', 16.831),
                                                                  ('PIC5', 15.207),
                                                                  ('PPC4', 1.752),
                                                                  ('PPC5', 1.079),
                                                                  ('foldover', 0),
                                                                  ('notfoldover', 1)]))

    def test_calculateConferencePareto_parallel(self):
        arrays = [oapackage.exampleArray(idx, 0) for idx in [45, 46, 47, 48]]

        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            presults, pareto = oaresearch.research_conference.calculateConferencePareto(
                arrays, N=None, k=None, verbose=0, number_parallel_jobs=2)
            stdout = mock_stdout.getvalue()
            self.assertEqual(stdout, 'Pareto: 2 optimal values, 2 objects\n')
        self.assertEqual(presults['nclasses'], 2)
        self.assertEqual(presults['pareto_indices'], (0, 3))
        self.assertEqual(presults['pareto_data'][0], OrderedDict([('ranksecondorder', 20),
                                                                  ('rankinteraction', 19),
                                                                  ('F4', (0, 1,
                                                                          4, 54, 11)),
                                                                  ('B4', 3.16),
                                                                  ('PEC4', 0.986),
                                                                  ('PEC5', 0.929),
                                                                  ('PIC4', 16.831),
                                                                  ('PIC5', 15.207),
                                                                  ('PPC4', 1.752),
                                                                  ('PPC5', 1.079),
                                                                  ('foldover', 0),
                                                                  ('notfoldover', 1)]))

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

        j4s = [oaresearch.research_conference.conferenceJ4(
            al) for al in arrays]
        self.assertEqual([np.sum(j4) for j4 in j4s], [52, 36, 36, 144])

        N = arrays[0].n_rows
        b4s = [np.sum(np.array(j4)**2) / N**2 for j4 in j4s]
        self.assertEqual(b4s, [3.16, 3.72, 3.4, 6.0])

        F4_values = [al.FvaluesConference(
            number_of_columns=4) for al in arrays]
        self.assertEqual(F4_values, [
                         (0, 1, 4, 54, 11), (0, 1, 8, 52, 9), (0, 2, 4, 51, 13), (0, 0, 27, 42, 1)])

    def test_conference_design_has_extensions(self):
        array = oapackage.exampleArray(42, 0)
        result = conference_design_has_extensions(array)
        self.assertEqual(result, True)

        array = oapackage.exampleArray(55, 0)
        result = conference_design_has_extensions(array)
        self.assertEqual(result, False)

        array = oapackage.exampleArray(55, 0)
        array = array.selectColumns([10, 11])
        result = conference_design_has_extensions(array)
        self.assertEqual(result, True)

        array = oapackage.exampleArray(55, 0)
        array = array.selectColumns([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11])
        result = conference_design_has_extensions(array)
        self.assertEqual(result, True)


if __name__ == '__main__':
    unittest.main()
