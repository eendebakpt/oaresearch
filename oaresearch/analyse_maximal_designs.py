#%%

# %% Load necessary packages
import os
import sys
import platform
import numpy as np
from importlib import reload
from os.path import join


import oapackage
import oapackage.graphtools
import oaresearch.research_conference
from oaresearch.research_conference import calculateConferencePareto, conferenceResultsFile, generateConferenceResults, \
    conferenceDesignsPage


# %%
    
def conference_design_extensions(array, verbose=0):
    """ Return list of ALL extensions """
    j1zero = 0
    conference_type = oapackage.conference_t(array.n_rows, array.n_columns, j1zero)

    zero_index = -1
    filterj2 = 1
    filterj3 = 0
    filter_symmetry = 1  # we can use symmetry reduction, since any the other filtering is not related to the symmetry of the design
    extensions = oapackage.generateSingleConferenceExtensions(
        array, conference_type, zero_index, verbose >= 2, filter_symmetry, filterj2, filterj3, filter_symmetry)
    extensions=[oapackage.hstack(array, extension) for extension in extensions]
    
    return list(extensions)

from functools import reduce
import operator

def flatten(data):
    if len(data)==0:
        return data
    return reduce(operator.concat, data)

def extent_list(extensions):
    return flatten([conference_design_extensions(array) for array in extensions] )
    
#%%
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-14-6-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-12-6-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-10-6-pareto.oa'
#afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-16-13-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-16-11-pareto.oa'
#afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-16-10-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-16-9-pareto.oa' # case!
#afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-20-16-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-14-6-pareto.oa' # case!
#afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-10-6-pareto.oa' 

lst = oapackage.readarrayfile(afile)
print(f'found {len(lst)} arrays')
#

from oaresearch.research_conference import conference_design_has_extensions

for idx, array in enumerate(lst):
    print(f'array {idx} in {afile}')
    
    has_extensions=conference_design_has_extensions(array)
    maximal=not has_extensions

    if maximal:
        continue
    print(f'  maximal {maximal}') 

    extensions0 = list(conference_design_extensions(array))
    extensions = extensions0


    N=array.n_rows
    
    for c in range(extensions[0].n_columns, N ):
        
        extensions2=extent_list(extensions)
        extensionsr=oapackage.selectConferenceIsomorpismClasses(extensions2, verbose=0)
        print(f'columns {c}->{c+1}: {len(extensions)}->{len(extensions2)}->{len(extensionsr)}')
        extensions=extensionsr
    if len(extensions)==0:
        print(f'found a case for array {idx} {array}!!!')
        break

#%%        
array.showarray()
    
#%%
from oaresearch.research_conference import conference_design_has_extensions

elist=[]
for idx, array in enumerate(lst):
    print(f'array {idx}')
    
    has_extensions=conference_design_has_extensions(array)
    maximal=not has_extensions

    if maximal:
        continue
    print(f'  maximal {maximal}') 

    extensions = list(conference_design_extensions(array))
    elist.append(extensions)

#%%
import functools
def smaller(a,b):
    return oapackage.compareLMC0(a,b)
rr=[]
for xx in elist:
    reduced = [oapackage.reduceConference(array) for array in xx]
    rr.append(sorted(reduced, key=functools.cmp_to_key(smaller)))
    

def equal_set(A,B):
    n=len(A)
    return np.all( [A[ii]==B[ii] for ii in range(n)])

ix=0
iy=1
for ii in range(len(rr[ix])):
    print(rr[ix][ii]==rr[iy][ii])

print(f'equal_set 0-1: {equal_set(rr[0], rr[1])}')


    