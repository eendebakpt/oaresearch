# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 10:12:51 2020

@author: eendebakpt
"""

#%% Load packages

import numpy as np
import logging
import itertools
import oapackage

#%% Define method to check irredundancy and format an array to a quantum state

def is_irredundant(array, k : int, verbose=0):
    """ Return True if the specified array is irredundant with respect to the number of columns specified
    
    For more details and the relation between irredundant orthogonal arrays and k–multipartite maximally
    entangled quantum states, see "Genuinely multipartite entangled states and orthogonal arrays", 2014, https://arxiv.org/pdf/1404.3586.pdf

    Args:
        array: Input array
        k: Check all tuples of number of columns - k columns of the array
        
    Returns:
        True is the array is k-irredundant

    """
    
    numpy_array = np.array(array)
    nruns = numpy_array.shape[0]
    ncolumns = numpy_array.shape[1]
    if k==ncolumns:
        return True
    for combination in (itertools.combinations(range(al.n_columns),al.n_columns-k)):
        logging.debug(f'check combination {combination}')

        #subarray = al.selectColumns(combination)
        #idx = oapackage.sortrows(subarray)

        np_subarray=numpy_array[:, combination]
        idx = np.lexsort( numpy_array[:,::-1].T, axis=0)
        
        for ii in range(nruns-1):
            
            if np.array_equal(np_subarray[idx[ii],: ],np_subarray[idx[ii+1],:]):
                if verbose:
                    print(f'array is not {k}-irredundant: columns {combination}, rows {idx[ii]}, {idx[ii+1]}')
                return False
    return True

def format_state(al ):
    """ Format array as a quantum state 
    
    Returns:
        String representation of the state
    """
    s=''
    pairs = []
    for row in np.array(al):
        vec=''.join([f'{d}' for d in row])
        pairs.append( f'|{vec}>' )
    s='+'.join(pairs)
    return s

    
al=oapackage.exampleArray(2,1)

al.showarray()
k=2
is_irredundant(al, k,1 )

#%%
#ll=oapackage.oahelper.runExtend(8, 4, t=2)

# ncols=2; nrows=2; strength=1 # generates the Bell state

ncols=3; nrows=4; strength=2 # 

#ncols=3; nrows=2; strength=1 # ! do not runm, this crashes due to limitations of oapackage, should generate GHZ for 3-qubits?

#ncols=4; nrows=8; strength=2 # generates Phi4 (equation (5) from https://arxiv.org/pdf/1404.3586.pdf)

#ncols=6; nrows=16; strength=2
#ncols=8; nrows=16; strength=2 # generates 2-irredundant states
#ncols=12; nrows=16; strength=2
#ncols=8; nrows=16; strength=3

ncols=8; nrows=16; strength=3

if strength==1:
    print('!! software was not written for strength 1, oapackage might crash')
class_name=f'OA({nrows}, 2^{ncols}, {strength})'
ll=oapackage.oahelper.runExtend(nrows, ncols, t=strength)

if len(ll)==0:
    print('no arrays generated')
for idx, al in enumerate(ll):
    r=False
    for  k in range(1, strength+1)[::-1]:
        r=is_irredundant(al, k=k, verbose=0 )
        if r:
            print(f'array {idx} in {class_name} is {k}-irredundant {r}')
            st=format_state(al)    
            print('state: ' + st)
            break
    if not r:
        print(f'array {idx} is not irredundant')
    


#%%

#ad=oapackage.arraydata_t(3, 9, 2, 4)

ll=oapackage.oahelper.runExtend(9, 3, t=2, l=3)
print(ll)
al=ll[0]
al.showarraycompact()
    

is_irredundant(al, k=2, verbose=1)
is_irredundant(al, k=1, verbose=1)