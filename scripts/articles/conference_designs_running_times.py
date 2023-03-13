"""

Example script for calculating conference designs

Pieter Eendebak <pieter.eendebak@gmail.com>
"""

# %% Load necessary packages

import platform
import sys
import time
from os.path import join

import numpy as np
import oapackage
import oapackage.graphtools

# %%


def generateConference(
    N, kmax=None, verbose=1, diagc=False, mode="nauty", nmax=None, selectmethod="random", tag="cdesign", outputdir=None
):
    """General function to compute conference matrices

    Arguments:
        N : integer
            number of rows in the array
        kmax : integer
            maximum number of columns to compute
        verbose : integer
            output level
        diagc : boolean
            the default value is False. If True, then only the diagonal
            matrices will be computed (e.g. all zeros are on the diagonal)

    """
    if kmax is None:
        kmax = N
    ctype = oapackage.conference_t(N, N, 0)

    if diagc:
        ctype.ctype = oapackage.conference_t.CONFERENCE_DIAGONAL
        tag += "-diagonal"
    if nmax is not None:
        tag += "-r"

    al = ctype.create_root()

    ll = oapackage.arraylist_t()
    ll.push_back(al)
    LL = [[]] * (kmax)
    LL[1] = ll
    if verbose:
        print(f"generateConference: start: {ctype}, mode {mode}")
    if outputdir is not None:
        _ = oapackage.writearrayfile(join(outputdir, "cdesign-%d-%d.oa" % (N, 2)), LL[1], oapackage.ATEXT, N, 2)

    for extcol in range(2, kmax):
        if verbose >= 2:
            print("generateConference: N %d, extcol %d: %d designs" % (N, extcol, len(LL[extcol - 1])))
            sys.stdout.flush()
        LL[extcol] = oapackage.extend_conference(LL[extcol - 1], ctype, verbose=verbose >= 2)

        if mode == "nauty":
            LL[extcol] = oapackage.selectConferenceIsomorpismClasses(LL[extcol], verbose >= 2)
        elif mode == "lmc0":

            LL[extcol] = oapackage.selectLMC0(LL[extcol], verbose >= 2, ctype)
        else:
            raise NotImplementedError(f"mode {mode}")

        LL[extcol] = oapackage.sortLMC0(LL[extcol])

        if nmax is not None:
            na = min(nmax, len(LL[extcol]))
            if na > 0:
                if selectmethod == "random":
                    idx = np.random.choice(len(LL[extcol]), na, replace=False)
                    LL[extcol] = [LL[extcol][i] for i in idx]
                elif selectmethod == "first":
                    LL[extcol] = [LL[extcol][i] for i in range(na)]
                else:
                    # mixed
                    raise Exception("not implemented")
        afmode = oapackage.ATEXT
        if len(LL[extcol]) > 1000:
            afmode = oapackage.ABINARY
        if outputdir is not None:
            _ = oapackage.writearrayfile(
                join(outputdir, "%s-%d-%d.oa" % (tag, N, extcol + 1)), LL[extcol], afmode, N, extcol + 1
            )

    ll = [len(l) for l in LL]
    if verbose:
        print(f"   N {N}: generated sequence: {ll}")
    return LL


# %% Measure generation time
outputdir = None
cases = range(4, 22, 2)

generation_times = {}
for NN in cases:
    t0 = time.time()
    LL = generateConference(N=NN, kmax=NN + 1, outputdir=outputdir, verbose=1, mode="lmc0")
    dt = time.time() - t0
    generation_times[f"cdesign{NN}"] = dt

for NN in cases:
    dt = generation_times[f"cdesign{NN}"]
    print(f"conference designs N={NN}: {dt:.2f} [s]")

print(f"platform processor: {platform.processor()}")
