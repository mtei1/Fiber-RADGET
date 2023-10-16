#! /usr/bin/env python3


import numpy as np
from math import sqrt, ceil
import matplotlib.pyplot as plt
from skimage.filters import (
    butterworth as bw,
    difference_of_gaussians as dog,
    gaussian as g,
)
from scipy.ndimage import uniform_filter1d as uf, median_filter as med
from multiprocessing import Pool
from tqdm import tqdm
from morph import adecompose,decompose, factor
from scipy.signal import butter, sosfilt, sosfreqz

from datetime import datetime as dt
import textwrap

from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt


def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype="low", analog=False)


def blf(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def bbf1(x):
    return bbf(*x)


def blf1(x):
    return blf(*x)


def med1(x):
    return med(*x)


def decompose1(x):
    return adecompose(*x)[-1]


def pool(a, f, args, n):
    odata = np.empty(a.shape)
    with Pool(processes=n) as P:
        pargs = [(a[:, i],) + args for i in range(a.shape[1])]
        res = list(tqdm(P.imap(f, pargs), total=len(pargs)))
    for i, r in enumerate(res):
        odata[:, i] = r
    return odata


def main(args):
    start = dt.now()
    for path in args.data:
        if path.endswith("npy"):
            a = np.load(path)
        else:
            a = np.genfromtxt(path, delimiter=",")
        if args.ecol is not None:
            a = a[:, : args.ecol + 1]
        if args.scol is not None:
            a = a[:, args.scol :]

    ok = ~np.isnan(a[0, :])
    n = len(ok)
    defs = []

    i0 = 0
    while i0 < n:
        while i0 < n and not ok[i0]:
            i0 = i0 + 1
        if i0 >= n or not ok[i0]:
            break
        i1 = i0 + 1
        while i1 < n and ok[i1]:
            i1 += 1
        defs.append((i0, i1))
        i0 = i1 + 1
    print(defs)

    vmin, vmax = np.nanmin(a), np.nanmax(a)
    vmin -= vmin % 100
    vmax -= vmax % 100
    vmax += 100
    levels = range(round(vmin), round(vmax), int(vmax - vmin) // 20)
    F = [("orig", a)]
    print("LEVELS",list(levels))
    print("Calculating filters")
    for s in args.lowpass:
        print("Lowpass",s)
        l, fs, o = [float(x) for x in s.split(",")]
        odata = np.empty(a.shape) + np.nan
        for (i0, i1) in defs:
            odata[:, i0:i1] = pool(a[:, i0:i1], blf1, (l, fs, o), args.nproc)
        F.append((f"BLF Low:{l},Ord:{o},Rate:{fs}", odata))

    for s in args.amorph:
        print("AMorph",s)
        hws = [x.split(":") for x in s.split(",")]
        ihws = (tuple((t, float(h), w) for (t, h, w) in hws),)
        odata = np.empty(a.shape) + np.nan
        for (i0, i1) in defs:
            print(i0,i1,ihws)
            odata[:, i0:i1] = pool(a[:, i0:i1], decompose1, ihws, args.nproc)
        F.append((f"AMORPH Size:{s}", odata))        
        
    for s in args.morph:
        print("Morph",s)
        hws = [x.split(":") for x in s.split(",")]
        ihws = (tuple((t, float(h), w) for (t, h, w) in hws),)
        #odata = np.empty(a.shape) + np.nan
        #for (i0, i1) in defs:
        #    odata[:, i0:i1] = pool(a[:, i0:i1], decompose1, ihws, args.nproc)
        odata = decompose(a, ihws[0])[-1]
        F.append((f"MORPH Size:{s}", odata))

    for s in args.median:
        print("Median",s)
        size = int(s)
        odata = np.empty(a.shape) + np.nan
        for (i0, i1) in defs:
            odata[:, i0:i1] = med(a[:, i0:i1], (size, 1))
        F.append((f"MED Size:{s}", odata))

    for s in args.gaussian:
        print("Gaussian",s)
        sig, w = [float(x) for x in s.split(",")]
        odata = np.empty(a.shape) + np.nan
        for (i0, i1) in defs:
            odata[:, i0:i1] = g(a[:, i0:i1], sigma=sig, truncate=w, channel_axis=1)
        F.append((f"G Sig:{sig},Wid:{w}", odata))

    for s in args.uniform:
        print("Uniform",s)
        odata = np.empty(a.shape) + np.nan
        for (i0, i1) in defs:
            odata[:, i0:i1] = uf(a[:, i0:i1], size=int(s), mode="nearest", axis=0)
        F.append((f"UF Size:{s}", odata))

    for path in args.read:
        print("Reading data from file",path)
        F.append((f"PATH: {path}", np.load(path)))
        
    G = [F[0]]
    for x in F[1:]:
        G.append(x)
        if args.xmed is not None:
            print("Calculating spatial median filter")
            odata = np.empty(x[1].shape) + np.nan
            for (i0, i1) in defs:
                odata[:, i0:i1] = med(x[1][:, i0:i1], (1, args.xmed))
            G.append((x[0] + f"_med{args.xmed}", odata))

    F = G
    if args.dump:
        for name, data in F:
            path = name.replace(" ", "_").replace(":", "")
            print("Dump data to path",path)
            np.save(path, data)

    n = len(F)
    r, c = factor(n)
    print("Start ploitting")
    if args.oneplot:
        for idx, f in enumerate(F):
            plt.plot(f[1], label=f[0])
        plt.legend()
    else:
        fig, ax = plt.subplots(ncols=c, nrows=r)
        for idx, f in enumerate(F):
            title, data = f
            i, j = divmod(idx, c)
            print(i, j, c, r)
            if r == 1:
                if args.scol == args.ecol and args.scol is not None:
                    ax[j].plot(data.reshape(-1), label=title)
                    ax[j].set_title(title)
                else:
                    ax[j].contourf(data, levels=levels)
                    if args.contour is not None and idx > 0:
                        ax[j].contour(data, levels=[args.contour], colors="r")
                        ax[j].set_title(title + f" (CONT:{args.contour})")
                    else:
                        ax[j].set_title(title)
            else:
                if args.scol == args.ecol and args.scol is not None:
                    ax[i, j].plot(data.reshape(-1), label=title)
                    ax[i, j].set_title(title)
                else:
                    ax[i, j].contourf(data, levels=levels)
                    if args.contour is not None and idx > 0:
                        ax[i, j].contour(data, levels=[args.contour], colors="r")
                        ax[i, j].set_title(title + f" (CONT:{args.contour})")
                    else:
                        ax[i, j].set_title(title)

    dur = round((dt.now() - start).total_seconds())
    doy = 1 + round((start - start.replace(day=1,month=1,hour=0,minute=0,second=0,microsecond=0)).total_seconds() // 86400)
    plt.gcf().set_size_inches(30, 18)
    plt.savefig(f"fig_{doy}_{start:%H%M%S}_{dur}.png", dpi=300)
    #plt.gcf().set_size_inches(10, 6)
    #plt.savefig("fig.png", dpi=200)
    plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    general = parser.add_argument_group("general")
    general.add_argument(
        "-n",
        "--nproc",
        type=int,
        help="number of threads",
        default=1
    )
    general.add_argument(
        "-1",
        "--oneplot",
        help="single 1d-plot (requires scol=ecol)",
        action="store_true",
    )
    general.add_argument(
        "-c",
        "--contour",
        type=int,
        help="draw contour"
    )
    general.add_argument(
        "-x",
        "--xmed",
        type=int,
        help="apply an median filter in the x-direction"
    )

    data = parser.add_argument_group("data")
    data.add_argument(
        "-e", "--ecol",
        type=int,
        help="last column to use"
    )
    data.add_argument(
        "-s",
        "--scol",
        type=int,
        help="first column to use"
    )
    data.add_argument(
        "-d",
        "--dump",
        help="dump generated data", action="store_true"
    )
    data.add_argument(
        "-r",
        "--read",
        nargs="*",
        help="read data for plotting",
        default = []
    )
    data.add_argument(
        "data",
        nargs="*",
        help="data path[s], at least one is needed"
    )    

    filters = parser.add_argument_group(
        "filters",
        description=textwrap.dedent(
            """\
For each filter type several filters can be defined by separating
the parameter-lists (one for each filter) by spaces.

Parameters for each filter are given as comma separated lists and
are described below.

Arguments are floating point numbers unless stated otherwise."""
        ),
    )
    filters.add_argument(
        "-l",
        "--lowpass",
        nargs="*",
        type=str,
        help="lowpass (butterworth) filter: low,samplingrate,order",
        default=[],
    )
    filters.add_argument(
        "-g",
        "--gaussian",
        nargs="*",
        type=str,
        help="gaussians: sigma,width(unit sigma)",
        default=[],
    )
    filters.add_argument(
        "-u",
        "--uniform",
        nargs="*",
        type=str,
        help="uniform filter: size (integer, number of samples)",
        default=[],
    )
    filters.add_argument(
        "-m",
        "--median",
        nargs="*",
        type=str,
        help="median filter: size (integer, number of samples)",
        default=[],
    )
    filters.add_argument(
        "-M",
        "--morph",
        nargs="*",
        type=str,
        help="morphological filter: t0:h0:s0,t1:h1:s1,... arguments are comma separated list of triples which in turn consist of colon separated values. Each triple consists of type of the structuring element f (one of sin,tri,step), the height and finally the size (integer, half-width of the element)",
        default=[],
    )
    filters.add_argument(
        "-A",
        "--amorph",
        nargs="*",
        type=str,
        help="alternative implementation of morphological filter: t0:h0:s0,t1:h1:s1,... arguments are comma separated list of triples which in turn consist of colon separated values. Each triple consists of type of the structuring element f (one of sin,tri,step), the height and finally the size (integer, half-width of the element) or two sizes separated by x (e.g. 360x15)",
        default=[],    
    )
    args = parser.parse_args()
    print("OPTIONS")
    for k,v in args.__dict__.items():
        print(f"\t{k:8} => {v}")
    main(args)
