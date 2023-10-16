import numpy as np
from matplotlib import pyplot as plt
from math import sin, ceil, sqrt

from scipy.ndimage import grey_dilation as dilation2d
from scipy.ndimage import grey_erosion as erosion2d
from scipy.ndimage import grey_opening as opening2d
from scipy.ndimage import grey_closing as closing2d


def factor(n):
    r, c = ceil(2 * sqrt(n / 6)), ceil(3 * sqrt(n / 6))
    if (r - 1) * c >= n:
        r -= 1
    if r * (c - 1) >= n:
        c -= 1
    return r, c


def g_sin(w, k):
    if len(w) == 1:
        l = w[0]
        a = np.array([sin(0.5 * np.pi * i / l) for i in range(0, 2 * l + 1)]).reshape(2*l+1,1)
    else:
        l,m = w
        a = np.zeros((2*l+1,2*m+1))
        for i in range(2*l+1):
            for j in range(2*m+1):
                v = sin(0.5 * np.pi * i/l) * sin(0.5 * np.pi * j/m)
                a[l-i,m-j] = a[l-i,m+j] = a[l+i,m-j] = a[l+i,m+j] = v
    return k * a

def g_tri(w, k):
    if len(w) == 1:
        l = w[0]
        a = np.array([i / l for i in range(0, l)] + [i / l for i in range(l, -1, -1)]).reshape(2*l+1,1)
    else:
        l,m = w
        a = np.zeros((2*l+1,2*m+1))
        for i in range(l+1):
            for j in range(m+1):
                a[l-i,m-j] = a[l-i,m+j] = a[l+i,m-j] = a[l+i,m+j] = (1-j/m)*(1-i/l)
    return k * a

def g_step(w, k):
    if len(w) == 1:
        l = w[0]
        a = np.ones((2*l+1,1))
    else:
        l,m = w
        a = np.ones((2*l+1,2*m+1))
    return k * a

def dilation(s, g):
    n, m = len(s), len(g)
    assert m < n
    a = np.zeros((n + m - 1, m), dtype=s.dtype)
    a += g
    for i in range(m):
        a[i : i + n, i] += s
        a[0:i, i] += s[0]
        a[i + n :, i] += s[-1]
    return np.max(a, axis=1)[m // 2 : m // 2 + n]


def erosion(s, g):
    n, m = len(s), len(g)
    assert m < n
    a = np.zeros((n + m - 1, m), dtype=s.dtype)
    a -= g
    for i in range(m):
        a[i : i + n, i] += s
        a[0:i, i] += s[0]
        a[i + n :, i] += s[-1]
    return np.min(a, axis=1)[m // 2 : m // 2 + n]


def closing(s, g):
    return closing2d(s, structure=g)


def opening(s, g):
    return opening2d(s, structure=g)


def aclosing(s, g):
    g = g.reshape(-1)
    return erosion(dilation(s,g),g)


def aopening(s, g):
    g = g.reshape(-1)
    return dilation(erosion(s,g),g)


def oc(s, g):
    return closing(opening(s, g), g)


def co(s, g):
    return opening(closing(s, g), g)

def aoc(s, g):
    return aclosing(aopening(s, g), g)


def aco(s, g):
    return aopening(aclosing(s, g), g)


def m(s, g):
    return (co(s, g) + oc(s, g)) / 2

def am(s, g):
    return (aco(s, g) + aoc(s, g)) / 2


def decompose(s, gs):
    fs = [s]
    for (t, h, w) in gs:
        if "x" in w:
            w = tuple([int(x) for x in w.split("x")])
        else:
            w = int(w),
        if t == "sin":
            g = g_sin(w, h)
        elif t == "tri":
            g = g_tri(w, h)
        else:
            g = g_step(w, h)
        fs.append(m(fs[-1], g))
    return fs


def adecompose(s, gs):
    fs = [s]
    for (t, h, w) in gs:
        if "x" in w:
            w = tuple([int(x) for x in w.split("x")])
        else:
            w = int(w),
        if t == "sin":
            g = g_sin(w, h)
        elif t == "tri":
            g = g_tri(w, h)
        else:
            g = g_step(w, h)
        fs.append(am(fs[-1], g))
    return fs


if __name__ == "__main__":
    from sys import argv

    path, idx = argv[1:3]
    if path.endswith("npy"):
        a = np.load(path)
    else:
        a = np.genfromtxt(path, delimiter=",")
    s = a[:, int(idx)]

    F = []
    for gpara in argv[3:]:
        gs = []
        for line in gpara.split(","):
            ws = line.split(":")
            gs.append((ws[0], float(ws[1]), int(ws[2])))
        fs = decompose(s, gs)
        F.append((gpara, fs))

    print(F)

    n = len(F) + 1
    r, c = factor(n)
    print(r, c)

    fig, ax = plt.subplots(ncols=c, nrows=r)
    for idx, f in enumerate(F):
        gpara, fs = f
        i, j = divmod(idx + 1, c)
        print(i, j, c, r)
        if r == 1:
            ax[j].set_title(gpara)
            ax[j].legend()
        else:
            ax[i, j].set_title(gpara)
            ax[i, j].legend()
        if r == 1:
            ax[0].plot(fs[-1], label=f[0])
            for idx, ff in enumerate(fs):
                ax[j].plot(ff, label=f"f_{idx}")
        else:
            ax[0, 0].plot(fs[-1], label=f[0])
            for idx, ff in enumerate(fs):
                ax[i, j].plot(ff, label=f"f_{idx}")
        if r == 1:
            ax[0].legend()
        else:
            ax[0, 0].legend()
    plt.savefig("plot.png")
    plt.show()
