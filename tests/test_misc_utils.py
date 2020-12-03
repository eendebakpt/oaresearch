# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

# %%
import numpy as np
import unittest

import oapackage
import oaresearch.research
from oaresearch.misc_utils import index_sorted_array


# %%


class TestMiscUtils(unittest.TestCase):

    def test_index_sorted_array(self):
        value = index_sorted_array( range(10),3)
        self.assertEqual(value, 3)

    def test_index_sorted_array_non_sorted_array(self):
        with self.assertRaises(ValueError):
            index_sorted_array( [1,4,1,1,2,2,2,2,3,3,5],4)
    
    def test_index_sorted_array_empty_array(self):
        with self.assertRaises(ValueError):
            index_sorted_array( [], 1)

    def test_index_sorted_array_no_valuesel(self):
        with self.assertRaises(ValueError):
            index_sorted_array( [1,2,3], 4)

if __name__ == '__main__':
    unittest.main()
