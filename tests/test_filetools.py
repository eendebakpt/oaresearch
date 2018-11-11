# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

import os
import unittest
import tempfile
import oapackage
import oaresearch.filetools

#%%


def create_oa_file(example_arrays=[1]):
    filename = tempfile.mktemp(suffix='.oa')
    oapackage.writearrayfile(filename, [oapackage.exampleArray(idx) for idx in example_arrays])
    return filename


class TestResearchoaFiletools(unittest.TestCase):

    def test_copyOAfile(self):
        source = create_oa_file([1, 1])
        narrays = oapackage.nArrays(source)
        self.assertEqual(narrays, 2)
        target = tempfile.mktemp(suffix='.oa')
        targetdir, target0 = os.path.split(target)
        resultfile = oaresearch.filetools.copyOAfile(source, targetdir, target0)
        self.assertEqual(resultfile, target0)

        target = tempfile.mktemp(suffix='.oa')
        targetdir, target0 = os.path.split(target)
        resultfile = oaresearch.filetools.copyOAfile(source, targetdir, target0, convert='B', zipfile=False)
        self.assertEqual(resultfile, target0)


if __name__ == '__main__':
    unittest.main()
