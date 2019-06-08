#%%

# %% Load necessary packages
import numpy as np
from importlib import reload
from os.path import join


import oapackage
import oapackage.graphtools
import oaresearch.research_conference



# %%
import oaresearch.research_conference
reload(oaresearch.research_conference)
from oaresearch.research_conference import maximal_extension_size, extend_conference_designs_full


maximum_number, designs = maximal_extension_size(oapackage.exampleArray(50,1))


#%%
import time
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-18-5-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-20-6-pareto.oa'
afile=r'C:/Users/eendebakpt/Dropbox/conference designs/designs/confpage_dc_extensions/conference/cdesign-24-11-pareto.oa'
afile=r'C:/Users/eendebakpt/Dropbox/conference designs/designs/confpage_dc_extensions/conference/cdesign-22-9-pareto.oa'

lst = oapackage.readarrayfile(afile)
array=lst[0]
t0=time.time()
maximum_number, designs = maximal_extension_size(array)
dt=time.time()-t0
print(f'dt {dt:.2f} [s]')

#%%
afile=r'C:/Users/eendebakpt/Dropbox/conference designs/designs/confpage_dc_extensions/conference/cdesign-22-7-pareto.oa'
lst = oapackage.readarrayfile(afile)

for array in lst:
    t0=time.time()
    maximum_number, designs = maximal_extension_size(array)
    dt=time.time()-t0
    print(f'dt {dt:.2f} [s]')


#%% Find specific cases
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

    if 0:    
        has_extensions=conference_design_has_extensions(array)
        maximal=not has_extensions
    
        if maximal:
            continue
        print(f'  maximal {maximal}') 

    maximum_number, designs = maximal_extension_size(array, verbose=0)

    print(f'array {idx}: maximum_number {maximum_number}, {len(designs)} designs')
    
    #extensions0 = list(conference_design_extensions(array))
    #extensions = extensions0


    N=array.n_rows
    if maximum_number==N:
        print(f'found a case for array {idx}: {array}!!!')
        break

#%% Full calculation
        
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-20-16-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-20-10-pareto.oa'
#afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-20-4-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-18-10-pareto.oa'
afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-18-5-pareto.oa'
#afile=r'C:/Users/eendebakpt/oatmp/confpage_dc_extensions/conference/cdesign-26-5-pareto.oa'
lst = oapackage.readarrayfile(afile)

for idx, array in enumerate(lst):
    print(f'array {idx} in {afile}')

    maximum_number, designs = maximal_extension_size(array, verbose=1)

    N=array.n_rows
    print(f'array {idx}: maximum_number {maximum_number}/{N}, {len(designs)} designs')
    

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


    