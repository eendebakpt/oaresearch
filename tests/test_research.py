# -*- coding: utf-8 -*-
"""
Pieter Eendebak <pieter.eendebak@gmail.com>

"""

#%%
import os
import numpy as np
import time
import itertools
import unittest

import oapackage
import oaresearch.research

#%%

def array2cpp(array, padding='   '):
    """ Convert array to C++ initialization code """    
    ss=padding + 'array_link array (%d, %d, 0);' % (array.n_rows, array.n_columns) +' \n' 
    ss+=padding +'int array_data_tmp[] = {%s};\n'  % (','.join(['%d'  % v for v in np.array(array).T.flatten()]))
    ss+=padding+'array.setarraydata (array_data_tmp, array.n_rows * array.n_columns);\n'
    return ss

class TestFunctions(unittest.TestCase):
    
    def test_array2cpp(self):
        al = oapackage.exampleArray(1)
        cpp=oaresearch.research.array2cpp(al)


if __name__=='__main__':
    unittest.main()
    