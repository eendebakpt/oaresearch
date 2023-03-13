#!/usr/bin/env python3
"""
Created on Wed May 18 12:52:15 2022

@author: eendebakpt
"""

# %%
import json
import os
import pathlib
import shutil
import subprocess
import tempfile
import time
from collections import namedtuple
from os.path import join

import matplotlib.pyplot as plt
import numpy as np
import oapackage
import rapidjson
from oapackage import splitFile
from oapackage.oahelper import arrayfile_generator
from qtt.utilities.tools import flatten, measure_time
from qtt.utilities.visualization import combine_legends
from termcolor import colored
from tqdm import tqdm

from oaresearch.pythondevelop import researchOA


def maxj_design(A):
    js = oapackage.jstruct_t(A.selectFirstColumns(5), 5)
    return js.maxJ()


def cprint(s):
    print(colored(s, color="cyan"))


def count_maxj(afile, outfile=None):
    cmd = f"oa_select_maxj -i {afile} -o {outfile} -f Z"
    r = os.system(cmd)
    if r != 0:
        raise Exception(f"count_maxj: return value {r}")
    nn = oapackage.nArrayFile(afile)
    nnm = oapackage.nArrayFile(outfile)
    return nn, nnm


def bprint(s):
    print(colored(s, color="cyan"))


def profile_expression(expression: str):
    """Profile an expression with cProfile and display the results using snakeviz"""
    import cProfile

    tmpdir = tempfile.mkdtemp()
    statsfile = os.path.join(tmpdir, "mystats")
    t0 = time.perf_counter()
    cProfile.run(expression, statsfile)
    dt = time.perf_counter() - t0
    print(f"profiling: {dt:.2f} [s]")
    # TODO: catch stdoutput of process
    r = subprocess.Popen(["snakeviz", statsfile])
    return r


def check_file(afile, tag=None):
    if isinstance(afile, str):
        aa = oapackage.readarrayfile(afile)
    else:
        aa = afile

    ngood = 0
    for ii, A in enumerate(aa):
        js = oapackage.jstruct_t(A, 5)
        mj = js.maxJ()
        if mj >= 32 + 1:
            ngood += 1

            js.Fval()

            A.Fvalues(5)

            # print(f'k {k}: {s1}.{s2}: ', ii, mj, J5)
    print(f"{tag}: {len(aa)} array(s), {ngood} have maxJ")
    return (len(aa), ngood)


# % Analyse starting point


xdir = "/media/eendebakpt/KONIJN/hopper/run64"

adfull = oapackage.readConfigFile(join(xdir, "oaconfig.txt"))
afile = join(xdir, "result-64.2-2-2-2-2.oa.gz")
aa = oapackage.readarrayfile(afile)
k = aa[0].n_columns
print(f"start : {len(aa)} array(s) with 5 columns")
for A in aa:
    J5 = A.Fvalues(5)
    js = oapackage.jstruct_t(A, 5)
    mj = js.maxJ()

    print(mj, J5)


def oa_result_file(ss, k, tag="extend"):
    idstr = adfull.reduceColumns(k).idstr()
    return oapackage.splitFile(ss) + f"-{tag}-" + idstr + ".oa"


#%%
k = 5
a7 = adfull.reduceColumns(k)
afile = join(xdir, f"eo-{a7.idstr()}.oa.gz")
# afile = join(xdir, f"result-{a7.idstr()}.oa.gz")
aa = oapackage.readarrayfile(afile)
print(f"start : {len(aa)} array(s) with {k} columns")

for idx, A in enumerate(aa):
    J5 = A.Fvalues(5)
    print(idx, J5)
    js = oapackage.jstruct_t(A, 5)
    mj = js.maxJ()
    J7 = A.Fvalues(7)
    js7 = oapackage.jstruct_t(A, 7)
    mj7 = js7.maxJ()

    if np.sum(np.abs(A.Jcharacteristics(5))) == 0:
        print(f"design {idx}: max(J5) {mj}  max(J7) {mj7} sum J7 {sum(J7)} gwlp {A.GWLP()}")

#%%
B = A.selectFirstColumns(5)


# %% Analyse start at k=7

maxj = []

maxblocks = set()

for s1 in range(780):

    afile = join(xdir, f"sp0-split-{s1}.oa.gz")
    aa = oapackage.readarrayfile(afile)
    k = aa[0].n_columns
    # print(f'k {k}: {s1}.{s2}: {len(aa)} array(s)')
    # J5 = A.Fvalues(5)
    # print(J5)
    for ii, A in enumerate(aa):
        js = oapackage.jstruct_t(A, 5)
        mj = js.maxJ()
        if mj >= 48 + 1:
            js.Fval()

            J5 = A.Fvalues(5)

            print(f"k {k}: {s1} design {ii}: maxj {mj}, J5 {J5} ")
            maxj.append(A)
            maxblocks.add(s1)

# maxj=maxj[-2:-1]
# maxj=maxj[:-3]

adata = adfull.reduceColumns(k)

xdir = "/mnt/data/tmp2/"
oapackage.writearrayfile("/mnt/data/tmp2/m64.oa", maxj)
# oapackage.writearrayfile(f'/mnt/data/tmp2/m64-{adata.idstr()}.oa', maxj)

# %% Branch factor
cdir = r"/mnt/data/tmp"

br = namedtuple("branch", ["columns", "n", "index"])


def generate_files(lvl, splits, adfull):

    adata = adfull.reduceColumns(splits[lvl].columns)
    maxk = splits[lvl + 1].columns

    w = oapackage.splitFile([b.index for b in splits[:lvl]])
    ifile_split = f"m64-{w}-{adata.idstr()}.oa"
    if w:
        tag_split = f"m64-{w}-sp{lvl}"
    else:
        tag_split = f"m64-sp{lvl}"
    w = oapackage.splitFile([b.index for b in splits[: lvl + 1]])

    infile_extend = f"m64-{w}.oa"  # -{adata.idstr()}.oa'
    outfile_extend = f"m64-{w}"
    return ifile_split, tag_split, infile_extend, outfile_extend, maxk


splits = [
    br(9, 100, 93),
    br(10, 30, 28),
    br(11, 50, 23),
    br(12, 2, 1),
    br(13, 2, 1),
    br(15, 2, 1),
    br(17, 2, 1),
    br(18, 1, 0),
]

splits = [br(13, 24 * 50, 3), br(14, 3, 2), br(15, 4, 3), br(18, 2, 1), br(22, 1, 0)]
lvl = 0


def oa_result_file(ss, k, tag="extend"):
    idstr = adfull.reduceColumns(k).idstr()
    return oapackage.splitFile(ss) + f"-{tag}-" + idstr + ".oa"


xdir = "/mnt/data/tmpx/sp0-split-549/sp1-split-6"
ifile0 = "sp0-split-549-sp1-split-6-extend-64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa.gz"

xdir = "/mnt/data/tmpx/sp0-split-549/sp1-split-190"
ifile0 = "sp0-split-549-sp1-split-190-extend-64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa.gz"
ifile_split, tag_split, infile_extend, outfile_extend, maxk = generate_files(0, splits, adfull)

print(f"narrays at start: {oapackage.nArrays(join(xdir, ifile0)):.3g}")

cprint("Stage")
if 0:
    ifile_splitx = ifile_split + ".gz" if ifile0.endswith("gz") else ifile_split
    cmd = f"cp {ifile0} {ifile_splitx}"
else:
    ifile_splitx = ifile_split
    cmd = f"oa_select_maxj -i {ifile0} -o {ifile_splitx} -f D"
print(cmd)

fid = open(os.path.join(cdir, "c.sh"), "w")
for lvl in range(0, len(splits) - 1):
    cprint(f"lvl: {lvl}")
    # b = splits[lvl]
    ifile_split, tag_split, infile_extend, outfile_extend, maxk = generate_files(lvl, splits, adfull)

    if lvl == 0:
        ifile_split = ifile_splitx
    cmd = f"oasplit -i {ifile_split} -n {splits[lvl].n} -f Z -o {tag_split}"
    cmd2 = (
        f"oaextendsingle -m 9 -l 2 -c ~/run64/oaconfig.txt  -f Z --maxk {maxk} -r {infile_extend} -o {outfile_extend}"
    )

    print(f"{cmd}; {cmd2}")

    fid.write(f"# lvl {lvl}\n{cmd}; {cmd2}\n")
fid.close()


# %% Analyse


br_result = namedtuple("br_result", ["ndesigns", "branch_factor", "k"])

xdir = "/mnt/data/tmpx/sp0-split-549/sp1-split-6"
bf = 1
prev_file = ifile0
rx = []
for lvl, b in enumerate(splits[:-2]):
    # print(b)
    bnext = splits[lvl + 1]
    bf = bf * splits[lvl].n

    ifile_split, tag_split, infile_extend, outfile_extend, maxk = generate_files(lvl, splits, adfull)
    next_file, _, _, _, _ = generate_files(lvl + 1, splits, adfull)

    print(ifile_split, next_file)
    nn0 = nn = oapackage.nArrayFile(join(xdir, ifile_split))
    if nn == -1:
        nn0 = nn = oapackage.nArrayFile(join(xdir, ifile_split) + ".gz")
    # break
    nn1 = oapackage.nArrayFile(join(xdir, infile_extend))
    nn = nn2 = oapackage.nArrayFile(join(xdir, next_file))
    cprint(f"{b.columns} -> {bnext.columns}: extend {nn0}->{nn1}->{nn}")
    # print(ifile)

    ndesigns = nn * bf
    print(f"  {bnext.columns}: branch factor {bf:.1f}, estimated: {ndesigns} = {ndesigns:.3g}")
    rx.append(br_result(ndesigns, bf, maxk))
c13 = {0: 0, 16: 5943257870, 32: 1371594144, 48: 11237518663, 64: 1433821720}

# c13=jd[13]

frac_start = rx[0][0] / c13[64]
dt = 11 * 60 + 3  # [s]

bfx = rx[1].branch_factor
r = rx[1]
tt = sum(jd[15].values())

maxj_frac = (jd[13][64]) / sum(jd[13].values())
print(ifile0)
print(splits)
print(f"estimated time: {bfx*dt/(24*3600):.2f} [days] (dt {dt:.1f} [s]")
print(f"estimated maxj designs at {r.k}: {r.ndesigns:.3g}/{tt:.2g} ")

# c15=

# %% Stored results
"""

m64--64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa m64-sp0-split-3-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa
13 -> 14: extend 5031940->4194->12708
  14: branch factor 1200.0, estimated: 15249600 = 1.52e+07
m64-sp0-split-3-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa m64-sp0-split-3-sp1-split-2-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa
14 -> 15: extend 12708->4236->4564
  15: branch factor 3600.0, estimated: 16430400 = 1.64e+07
m64-sp0-split-3-sp1-split-2-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa m64-sp0-split-3-sp1-split-2-sp2-split-3-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa
15 -> 18: extend 4564->1141->0
  18: branch factor 14400.0, estimated: 0 = 0
sp0-split-549-sp1-split-6-extend-64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa.gz
[branch(columns=13, n=1200, index=3), branch(columns=14, n=3, index=2), branch(
    columns=15, n=4, index=3), branch(columns=18, n=2, index=1), branch(columns=22, n=1, index=0)]
estimated time: 38.54 [days]
estimated maxj designs at 15: 1.64e+07/3.8e+10


m64--64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa m64-sp0-split-3-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa
13 -> 14: extend 5031940->4194->12708
  14: branch factor 1200.0, estimated: 15249600 = 1.52e+07
m64-sp0-split-3-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa m64-sp0-split-3-sp1-split-2-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa
14 -> 15: extend 12708->4236->4564
  15: branch factor 3600.0, estimated: 16430400 = 1.64e+07
m64-sp0-split-3-sp1-split-2-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa m64-sp0-split-3-sp1-split-2-sp2-split-3-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa
15 -> 18: extend 4564->1141->0
  18: branch factor 14400.0, estimated: 0 = 0
sp0-split-549-sp1-split-200-extend-64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa.gz
[branch(columns=13, n=1200, index=3), branch(columns=14, n=3, index=2), branch(
    columns=15, n=4, index=3), branch(columns=18, n=2, index=1), branch(columns=22, n=1, index=0)]
estimated time: 27.62 [days] (dt 663.0 [s]
estimated maxj designs at 15: 1.64e+07/3.8e+10


"""

# %% Extract all maxJ files (only once!)


aa = oapackage.readarrayfile(afile)

maxj_list = list(filter(maxj, aa))

# %%
for jj, A in enumerate(aa):
    print(jj)
    maxj(A)

# %%


def maxj_(A):
    js = oapackage.jstruct_t(A.selectFirstColumns(5), 5)
    js.maxJ()


class MaxJ:
    """Helper class to quickly count number of maxj designs"""

    def __init__(self, A):

        self.js = oapackage.jstruct_t(A.selectFirstColumns(5), 5)
        self.previous_array = A.selectFirstColumns(5) - 10
        self.previous_result = None

    def maxj(self, A):
        A5 = A.selectFirstColumns(5)

        if self.previous_array == A5:
            return self.previous_result

        self.previous_array = A5
        js = self.js  # oapackage.jstruct_t(A5, 5)
        js.calcj5(A5)
        self.previous_result = js.calculateF()[0]

        return self.previous_result


# %% Only once!


def deep_set(d, ss, value):
    if len(ss) == 1:
        d[ss[0]] = value
    else:
        deep_set(d[ss[0]], ss[1:], value)


def calculate_data(afile, outfile, tag, maxj=64):
    jsonfile = outfile.replace(".oa", ".json")
    if os.path.exists(jsonfile):
        jdata = json.load(open(jsonfile))
        nn = sum(jdata.values())

        nnm = jdata.get(str(maxj), 0)
        bprint(f"block {tag}: found {nnm}/{nn} maxj designs! (cache)")
    else:
        nn, nnm = count_maxj(afile, outfile=outfile)
        bprint(f"block {tag}: found {nnm}/{nn} maxj designs!")
    return nn, nnm


splits = [780, 248, 18]
# splits=[100, 4, 18]

adfull = oapackage.readConfigFile("/home/eendebakpt/run64/oaconfig.txt")
tmpdir0 = "/mnt/data/tmpx"
tmpdir0 = pathlib.Path(tmpdir0)
tmpdir0.mkdir(exist_ok=True)

odir = pathlib.Path("/mnt/data/maxjout")
odir.mkdir(exist_ok=True)

k = 13
adata = adfull.reduceColumns(k)

xdir = "/media/eendebakpt/KONIJN/hopper/run64"

l1_blocks = {545, 546}  # , 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558}:
l1_blocks = range(0, 780)
l1_blocks = [246]
# double special: block 567


odir = pathlib.Path("/mnt/data/maxjout")
odir_str = str(odir)

# %%


def unzip(s1, xdir, target_dir, run=True):
    """Unzip a stored file"""
    zfile = join(xdir, f"sp0-split-{s1}", f"splitted-{s1}.zip")

    sdir = os.path.join(target_dir, f"sp0-split-{s1}")
    oapackage.mkdirc(sdir)
    cmd = f"cd {target_dir}; unzip {zfile} -d {sdir}"
    print(cmd)
    if run:
        r = os.system(cmd)
    return r


# %%
counts = {}
counts_maxj = {}
mj_calculator = None


for s1 in l1_blocks:
    print(f"block {s1}")
    # copy and unzip data file
    zfile = join(xdir, f"sp0-split-{s1}", f"splitted-{s1}.zip")
    tmpdir = tmpdir0
    try:
        shutil.rmtree(tmpdir)
    except:
        pass
    tmpdir.mkdir(exist_ok=True)

    if s1 == 567:
        # not sure why 567 is not zipped...
        five_dir = join(xdir, f"sp0-split-{s1}")
        cmd = f"rsync -avr {five_dir}/ {tmpdir}/"
        print(cmd)
        # r = os.system(cmd)

    else:
        cmd = f"cd {tmpdir}; unzip {zfile} -d {tmpdir}"
        print(cmd)
        r = os.system(cmd)

    # loop over second level
    for s2 in range(splits[1]):

        idstr = adata.idstr()
        # print(afile)
        sdir = join(tmpdir, splitFile([0, s2])[12:])

        counts.setdefault(s1, {})
        counts_maxj.setdefault(s1, {})

        zfile0 = f"splitted-{oapackage.splitTag([s1, s2])}.zip"
        zfile = join(sdir, zfile0)
        if os.path.exists(zfile):
            print(f"block {s1}.{s2}: special case")
            cmd = f"cd {sdir}; unzip -o {zfile} -d {sdir}"
            print(cmd)
            r = os.system(cmd)
            assert r == 0

            counts[s1][s2] = {}
            counts_maxj[s1][s2] = {}
            for s3 in range(splits[2]):
                sdir2 = join(sdir, f"sp2-split-{s3}")

                ss = [s1, s2, s3]
                tag = "-".join(map(str, ss))
                # print(f'block {tag}: {nn} arrays')

                afile0 = oa_result_file(ss, k)
                afile = join(sdir2, afile0)

                outfile = join(odir, f"{tag}-k{k}-maxj64.oa")
                nn, nnm = calculate_data(afile, outfile, tag=".".join(map(str, ss)))
                deep_set(counts, ss, nn)
                deep_set(counts_maxj, ss, nnm)

        else:
            ss = [s1, s2]
            if s1 == 567:
                jfile_567 = os.path.join(five_dir, splitFile([0, s2])[12:], researchOA.numbersFile(ss, tag="jstats"))
                jdata = oapackage.readStatisticsFile(jfile_567, verbose=0)

                jd = {j: jdata.getCount(k, j) for j in [0, 16, 32, 48, 64]}

                jfile = "-".join([str(e) for e in ss]) + f"-k{k}-maxj64.json"
                json.dump(jd, open(join(odir_str, jfile), "w"))

                nn = sum(jd.values())
                nnm = jd.get(64, 0)
                print(f"block {ss}: {nn}, {nnm}")
            else:

                tag = "-".join(map(str, ss))
                afile0 = oa_result_file(ss, k)
                afile = join(sdir, afile0)

                outfile = join(odir, f"{tag}-k{k}-maxj64.oa")
                nn, nnm = calculate_data(afile, outfile, tag=".".join(map(str, ss)))
                # print(f'block {tag}: {nn} arrays')

            deep_set(counts, ss, nn)
            deep_set(counts_maxj, ss, nnm)
        # detect special....


# %%
json.dump((counts, counts_maxj), open(join(odir, f"counts-{k}.json"), "w"))

# %% Load pre-calculated data

(counts, counts_maxj) = json.load(open(join(odir, f"counts-{k}.json")))


# %%


def dict_sum(dd):
    all_keys = list(set(flatten([list(d.keys()) for d in dd])))

    return {k: sum(w.get(k, 0) for w in dd) for k in all_keys}


def getnumber(entry):
    """Get total number in subblock"""
    if isinstance(entry, dict):
        return sum(getnumber(v) for v in entry.values())
    else:
        return entry


def get_full_numbers(entry, k, start_block=[], subblocks=[]):
    """Get total number in subblock"""
    for b in start_block:
        entry = entry[str(b)]
    if isinstance(entry, dict):
        blocks = start_block + subblocks
        sub_results = [get_full_numbers(v, k, start_block=[], subblocks=blocks + [key]) for key, v in entry.items()]
        return dict_sum(sub_results)
    else:
        jfile = "-".join([str(e) for e in start_block + subblocks]) + f"-k{k}-maxj64.json"
        entry = rapidjson.load(open(join(odir_str, jfile), "rb"))
        # print(f' {start_block+subblocks}: {entry}')
        return entry


with measure_time():
    r = get_full_numbers(counts, k, [])
    print(r)

# %%

jdata = oapackage.readStatisticsFile((join(xdir, "jstats-base.txt")), verbose=0)

jd = {k: {j: jdata.getCount(k, j) for j in [0, 16, 32, 48, 64]} for k in range(jdata.maxCols() + 1)}
print(jd)

plt.figure(1)
plt.clf()

kk = range(8, 22)
for j in [16, 32, 48, 64]:
    xx = np.array([jd[k][j] if jd[k][j] >= 0 else np.NaN for k in kk])
    xx = np.maximum(xx, 1e-10)
    # print(xx)
    label = f"max($J_5$)={j}"
    if j == 64:
        idx = kk.index(13)
        h = plt.plot(kk[: idx + 1], xx[: idx + 1], ".-", label=label)[0]
        plt.plot(kk[idx:], xx[idx:], ".--", alpha=0.72, color=h.get_color())
    else:
        plt.plot(kk, xx, ".-", label=label)
plt.yscale("symlog")
plt.legend()
plt.xlabel("Number of columns")
plt.ylabel("Number of designs")
plt.title(f"Number of designs in series ${adfull.latexstr()}$")

# %%

jfile = "-".join([str(e) for e in start_block + subblocks]) + f"-k{k}-maxj64.json"


# %%
cc = []
for s1 in range(splits[0]):
    if str(s1) in counts:
        for s2 in range(splits[1]):
            num = getnumber(counts[str(s1)][str(s2)])
            if num == 0:
                raise
            cc.append(num)
ccm = []
tt = []
for s1 in range(splits[0]):
    if str(s1) in counts:
        for s2 in range(splits[1]):
            tt.append((s1, s2))
            ccm.append(getnumber(counts_maxj[str(s1)][str(s2)]))
v1 = getnumber(counts)
v2 = getnumber(counts_maxj)
print(v1, v2)

sidx = np.argsort(cc)
cc = [cc[idx] for idx in sidx]
ccm = [ccm[idx] for idx in sidx]

bins = np.logspace(np.log10(max(1, np.floor(np.min(cc)))), np.log10(np.max(cc)), 40)

bin_indices = np.digitize(cc, bins)

# %%
nblocks = [
    len(cc[np.searchsorted(bin_indices, j, side="left") : np.searchsorted(bin_indices, j, side="right")])
    for j in range(1, len(bins))
]
ndesigns = np.array(
    [
        sum(cc[np.searchsorted(bin_indices, j, side="left") : np.searchsorted(bin_indices, j, side="right")])
        for j in range(1, len(bins))
    ]
)
ndesigns64 = np.array(
    [
        sum(ccm[np.searchsorted(bin_indices, j, side="left") : np.searchsorted(bin_indices, j, side="right")])
        for j in range(1, len(bins))
    ]
)


# %%
# v=json.load(open(r'/home/eendebakpt/projects/oapackage/build/x.json'))
# print(v)


bin_centres = (bins[:-1] + bins[1:]) / 2
plt.figure(30)
plt.clf()
ax = plt.gca()

vv, bins, _ = plt.hist(cc, bins, label="Designs")

# plt.plot(cc, np.array(ccm)/np.array(cc), '.b')
# plt.hist(ccm, bins, color='r', rwidth=.68)
plt.xlabel("Number of designs in block")
plt.ylabel("Number of blocks")
# plt.semilogx()
plt.loglog()
ax2 = ax.twinx()
ax2.plot(bin_centres, np.array(ndesigns64 / ndesigns), ".-r", label="Fraction")
plt.ylabel("Fraction of designs with max(J5)=64")
plt.title(f"Distribution of designs over {len(cc)} blocks")
# plt.suptitle(f'(total number of blocks: {780*248})')
# plt.semilogx()
plt.legend()
combine_legends([ax, ax2])
g = plt.grid(which="major", axis="y")

# %%
plt.figure(10)
plt.clf()
vv, bins, _ = plt.hist(cc, bins, label="Designs")
plt.hist(ccm, bins, color="r", rwidth=0.68)
plt.xlabel("Number of designs")
plt.ylabel("Frequency")
plt.semilogx()
plt.title(f"Distribution of designs over {len(cc)} blocks with a max J design")
plt.suptitle(f"(total number of blocks: {780*248})")
# plt.semilogx()
# plt.legend()

plt.figure(20)
plt.clf()
plt.plot(cc, np.array(ccm) / np.array(cc), ".b")
# plt.hist(ccm, bins, color='r', rwidth=.68)
plt.xlabel("Number of designs")
plt.ylabel("Fraction of designs with max(J5)=64")
plt.semilogx()
plt.title(f"Distribution of designs over {len(cc)} blocks")
plt.suptitle(f"(total number of blocks: {780*248})")
# plt.semilogx()
# plt.legend()

# %%

vals = np.random.random(1e8)
vals.sort()

nbins = 100
bins = np.linspace(0, 1, nbins + 1)
ind = np.digitize(vals, bins)

results = [
    func(vals[np.searchsorted(ind, j, side="left") : np.searchsorted(ind, j, side="right")]) for j in range(1, nbins)
]


# %% Maxj designs at k columns

k = 13
for s1 in range(545, 559):
    r = get_full_numbers(counts, k, [s1])
    print(f"block {s1}: k {k}: {r}")


# So 546 and 549 are the largest blocks. We only need to consider one of them

for s2 in range(248):
    s = [549, s2]
    r = get_full_numbers(counts, k, s)
    r = {k: r[k] for k in sorted(r)}
    print(f"block {s}: k {k}: {r}")

# %%
unzip(549, xdir, tmpdir0)


# %%

s = [549, 6]
x = os.path.join(tmpdir0, oapackage.splitDir(s))
print(x)

afile = oa_result_file(s, 13)
ifile = os.path.join(x, afile + ".gz")
print(ifile)

shutil.copyfile(ifile, os.path.join(tmpdir0, "m64-64.2-2-2-2-2-2-2-2-2-2-2-2-2-2.oa"))


# %% TODO:
# Check block 567

# %% Histogram

values = [list(v.values()) for v in counts.values()]
values = np.array(values).flatten()

plt.figure(10)
plt.clf()
plt.hist(values, bins=40)


# %%


s1 = 549
s1 = 545


for s2 in range(248):

    xdir = f"/mnt/data/tmp/{s1}/sp1-split-{s2}"

    k = 12
    kk = "-".join(["2" for j in range(k)])
    afile = join(xdir, f"sp0-split-{s1}-sp1-split-{s2}-pareto-64.{kk}.oa.gz")
    aa = oapackage.readarrayfile(afile)
    tag = f"k {k}: {s1}.{s2}"
    # _,n12=check_file(aa, tag=tag)

    k = 14
    kk = "-".join(["2" for j in range(k)])
    afile = join(xdir, f"sp0-split-{s1}-sp1-split-{s2}-pareto-64.{kk}.oa.gz")
    aa = oapackage.readarrayfile(afile)
    tag = f"k {k}: {s1}.{s2}"
    _, n12 = check_file(aa, tag=tag)

    # for ii, A in enumerate(aa):
    #     js=oapackage.jstruct_t(A, 5)
    #     mj=js.maxJ()
    #     if mj>=32+1:
    #         js.Fval()

    #         J5 = A.Fvalues(5)

    #         print(f'k {k}: {s1}.{s2}: ', ii, mj, J5)

# %%


split = [545, 242]
split = [545, 87]
oapackage.splitFile(split)

k = 14
kk = "-".join(["2" for j in range(k)])


for snext in range(18):

    xx = oapackage.splitDir(split + [snext])
    xx = "/".join(xx.split("/")[1:])
    xdir = join("/mnt/data/tmp/", f"{split[0]}", xx)

    afile = join(xdir, oapackage.splitFile(split + [snext]) + f"-extend-64.{kk}.oa.gz")

    aa = oapackage.readarrayfile(afile)
    # print(f'k {k}: {s1}.{s2}: {len(aa)} array(s)')

    tag = f"k {k}: {oapackage.splitTag(split+[snext])}"
    check_file(aa, tag=tag)

# %%

s1 = 548
s2 = 91

k = 12

for k in range(10, 16):
    kk = "-".join(["2" for j in range(k)])

    xdir = f"/mnt/data/tmp/{s1}/sp1-split-{s2}"
    afile = xdir + "/" + f"sp0-split-{s1}-sp1-split-{s2}-extend-64.{kk}.oa.gz"

    aa = oapackage.readarrayfile(afile)

    ngood = 0
    for ii, A in enumerate(aa):
        js = oapackage.jstruct_t(A, 5)
        mj = js.maxJ()
        if mj >= 32 + 1:
            ngood += 1

            js.Fval()

            J5 = A.Fvalues(5)

            # print(f'k {k}: {s1}.{s2}: ', ii, mj, J5)
    print(f"k {k}: {s1}.{s2}: {len(aa)} array(s), {ngood} have maxJ")

# %% Time of L5 check
s = [567, 90]
s = [556]


oaextend = oapackage.OAextend()


def calculate_lmc_speed(s, niter=5, select_maxj=False):
    results = {}
    xdir = join("/home/eendebakpt/run64", oapackage.splitDir(s))
    for kk in range(9, 21):
        adata = adfull.reduceColumns(kk)
        afile = oapackage.splitFile(s) + "-pareto-" + adata.idstr() + ".oa"

        aa = oapackage.readarrayfile(join(xdir, afile))

        # print(f'{kk}: {len(aa)}')

        na = len(aa)
        aamax = list(filter(lambda x: maxj_design(x) != 64, aa))
        # print(f' maxj ratio: {len(aamax)}/{len(aa)}')
        if select_maxj:
            aa = aamax
        if len(aa) > 0:

            with measure_time(None) as mt:
                reduction = oapackage.LMCreduction_t(adata)
                for ij in range(niter):
                    for A in aa:
                        # oapackage.LMC_LESS = 0
                        reduction.reset()
                        oapackage.LMCcheckj5(A, adata, reduction, oaextend)
            dtn = (mt.dt / len(aa)) / niter
            print(f"k {kk}: ndesigns {len(aamax)}/{na}: {1/dtn:.4f} checks/s")
            results[kk] = {"n": len(aa), "dt": dtn, "cps": 1 / dtn}
    return results


results = calculate_lmc_speed(s)
print(results)

# %%
ad = []
# for s1 in range(20, 780, 10):
for s1 in range(545, 560, 1):
    # for s2 in [0, 20, 40, 80]:
    s = [s1]
    bprint(f"determine speed on {s}")
    results = calculate_lmc_speed(s, select_maxj=True)
    ad.append(results)

# %%
kk = range(8, 21)
speed = [None] * len(kk)
for ii, k in enumerate(kk):
    speed[ii] = np.nanmean([r.get(k, {"cps": np.NaN})["cps"] for r in ad])

plt.figure(20)
plt.clf()
plt.plot(kk, speed, ".-")
plt.xlabel("Number of columns")
plt.ylabel("Checks / second")
plt.suptitle("Performance of L5 check", fontsize=18)
# plt.title('(measured on selection of designs in L5 format)')
plt.title("(measured on selection of max($J_5$)=64 designs in L5 format)")

# %%

afile = join(xdir, r"test-64.2-2-2-2-2-2-2-2-2-2-2-2-2.oa")
check_file(afile)


# %%
arraylist = maxj
al = maxj[0]
adfull = oapackage.readConfigFile("/home/eendebakpt/run64/oaconfig.txt")
discardJ5 = 20
dextend = oapackage.depth_extend_t(adfull, 10, discardJ5)

j5structure = 1
oaextendx = oapackage.OAextend()
oaextendx.setAlgorithm(9, adfull)
oaextendx.j5structure = j5structure

oaextendx.checkarrays = 0
oaextendx.use_row_symmetry = 0
oaextendx.extendarraymode = oapackage.OAextend.extendarray_mode_t_APPENDFULL
oaextendx.init_column_previous = oapackage.INITCOLUMN_J5
oaextendx.nLMC = 100000
oaextendx.info()

# %%
maxarrays = 2e6
dextend.oaextend = oaextendx
dextend.logtime = 10
dextend.setNarraysMax(maxarrays)

initcolumn = adfull.strength + 1
initcolumn = arraylist[0].n_columns

#        dextend.arraywriter = new arraywriter_t ();
#        dextend.arraywriter->initArrayFiles (*dextend.ad, initcolumn + 1, outputprefix, mode);
#        dextend.arraywriter->writearrays = writedeptharrays;

dextend.counter = oapackage.counter_t(adfull.N)
dextend.loglevelcol = 7
# %%
verbose = 2
# depth_extensions_storage_t *ds = 0;
ai = 0

dextendloop = oapackage.depth_extend_t(dextend)
dextendloop.setposition(al.n_columns, ai, len(arraylist), 0, 0)

oapackage.depth_extend_array(al, dextendloop, adfull, verbose, None, ai)

# %%


profile_expression("go()")

# %%

# %% Legacy helper files


def count_maxj2(afile, outfile=None):
    """Count number of maxj elements in a file and store them in result file"""
    nn = oapackage.nArrayFile(afile)
    afout = None
    mj_calculator = None

    nnm = 0
    for A in tqdm(arrayfile_generator(afile), total=nn, leave=False):
        if mj_calculator is None:
            mj_calculator = MaxJ(A)
        if mj_calculator.maxj(A):
            if afout is None and outfile is not None:
                afout = oapackage.arrayfile_t(outfile, 64, A.n_columns, -1, oapackage.ABINARY_DIFFZERO, 1)
            if afout is not None:
                afout.append_array(A)
            nnm = nnm + 1
    if afout:
        afout.closefile()
    return nn, nnm
