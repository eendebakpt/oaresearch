"""
Created on Sat Sep 28 16:43:36 2013

@author: eendebakpt
"""

# %% Exammple of delete-one-factor article appendix
import oapackage

al = oapackage.exampleArray(4)
al = oapackage.reduceDOPform(al)
al.showarray()
print("GWLP %s" % str(al.GWLP()))
for ii in range(0, al.n_columns):
    bl = al.deleteColumn(ii)
    print("dof %d: GWLP %s" % (ii, str(bl.GWLP())))

dopgwp = oapackage.projectionGWLPs(al)
sg = oapackage.symmetry_group(dopgwp, False)
sg.show(1)

al = oapackage.exampleArray(5)
al.showarray()
print("GWLP %s" % str(al.GWLP()))

aldop = oapackage.reduceDOPform(al)
aldop.showarray()

# import oahelper
for ii in range(0, al.n_columns):
    bl = al.deleteColumn(ii)
    print("delete column %d: GWLP %s" % (ii, str(bl.GWLP())))

aldop = oapackage.reduceDOPform(al, 1)
