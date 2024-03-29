"""
Created on Fri Feb  8 22:37:26 2019

@author: eendebakpt
"""


# %%
import numpy
import numpy as np
import oapackage

array = oapackage.exampleArray(10, 1).selectFirstColumns(1)
array = oapackage.array_link(np.array(array)[::3, :])
# array=oapackage.exampleArray(0,1).selectFirstColumns(1);
al = array
arrayclass = oapackage.arraylink2arraydata(al, 0, 0)
s = arrayclass.factor_levels()

r = oapackage.array2eigenModelMatrixMixed(array, verbose=3)
N = al.n_rows
k = al.n_columns
print(r)

# %%
myprintf = print
df = np.array(s) - 1

myprintf("df ")
print(df)

AA = np.array(al)

# %%
mesize = numpy.sum(df)

print("main effects: size %dx%d\n" % (N, mesize))

# %%
npx = oapackage.numberModelParams(al)
main_effects = numpy.zeros((N, mesize))

# %%


def helmert_contrasts(number_of_levels, verbose=0):
    """Calculate Helmert contrasts for a given number of levels for a number_of_levels

    Args:
        number_of_levels: number of levels in the number_of_levels
    Returns:
        array: array with calculated Helmert contrasts
    """
    N = number_of_levels
    meoffset = 0
    md = number_of_levels - 1

    main_effects = numpy.zeros((number_of_levels, md))
    Z = np.zeros((number_of_levels, md + 1))

    for value in range(number_of_levels):
        for ii in range(0, md + 1):
            Z[value, ii] = value > ii - 1

        # make Helmert contrasts (these are automatically orthogonal)
        Z[value, 0] = 1
        if value > 0:
            Z[value, value] = value

        for q in range(1, value):
            Z[value, q] = 0
        for q in range(value + 1, md + 1):
            Z[value, q] = -1

    if verbose:
        print("helmert_contrasts: %d\n" % (number_of_levels))
        print("Z (initial creation)")
        print(Z)

    # normalize the contrasts
    for ii in range(0, md):
        normalization = Z[:, (ii + 1)].T.dot(Z[:, ii + 1])
        if verbose:
            print("helmert_contrasts: normalize number_of_levels tmp: %s " % (normalization,))
        main_effects[:, meoffset + ii : (meoffset + ii + 1)] = (
            np.sqrt(N) * Z[:, (ii + 1) : (ii + 2)] / np.sqrt(normalization)
        )

    return main_effects


value = 0
number_of_levels = arrayclass.getfactorlevel(0)
hc = helmert_contrasts(arrayclass.getfactorlevel(0), verbose=0)
print(hc)
# helmert_contrast(2, arrayclass.getnumber_of_levelslevel(0))


# %%
import oapackage

array = oapackage.exampleArray(0, 1)
array.showarray()
M = oapackage.array2modelmatrix(array, "i")
print(M)

# %%
import numpy as np

np.set_printoptions(precision=3)

array = oapackage.exampleArray(10, 1)

s = oapackage.array2modelmatrix_sizes(array)
print(s)

M = oapackage.array2modelmatrix(array, "i")
print(np.array(array)[0::, : s[1]])
print(M[0::, : s[1]])

# %%
# M=oapackage.array2modelmatrix(array, 'q')
# print(M[0:10,:])
