#!/usr/bin/env python3
"""count_basis2.py: TROVE contracted product-basis size for CH3OH vs ranges/polyad weights.
Classes: CO(1), OH(1), CH3(3), COH(1), OCH(3), S(2), tors(1). Class polyad = int(w*sum v) (int per class)."""
import itertools
from collections import Counter

def cls(nm, vmax, w):
    out = Counter()
    for v in itertools.product(range(vmax + 1), repeat=nm):
        if sum(v) <= vmax:
            out[int(w * sum(v) + 1e-9)] += 1
    return out

def tors(nf, w):
    out = Counter()
    for n in range(nf):
        out[int(w * n + 1e-9)] += 1
    return out

def count(r, w, nt, wt, pmax):
    cs = [cls(1, r['CO'], w['CO']), cls(1, r['OH'], w['OH']), cls(3, r['CH'], w['CH']),
          cls(1, r['COH'], w['COH']), cls(3, r['OCH'], w['OCH']), cls(2, r['S'], w['S']), tors(nt, wt)]
    dist = Counter({0: 1})
    for c in cs:
        nd = Counter()
        for p0, n0 in dist.items():
            for p, n in c.items():
                if p0 + p <= pmax:
                    nd[p0 + p] += n0 * n
        dist = nd
    return sum(dist.values())

w = dict(CO=1, OH=2, CH=2, COH=1, OCH=1, S=1)
if __name__ == "__main__":
    r4 = dict(CO=4, OH=4, CH=4, COH=4, OCH=4, S=4)
    r46 = dict(CO=4, OH=4, CH=4, COH=4, OCH=4, S=6)
    for name, r in (("v<=4 (S<=4)", r4), ("v<=4 (S<=6)", r46)):
        for nt, wt in ((9, 0.1), (9, 0.34), (12, 0.34), (15, 0.34)):
            for pmax in (6, 7, 8, 9, 10):
                n = count(r, w, nt, wt, pmax)
                print(f"{name:12s} tors {nt:2d} (w={wt}) Pmax={pmax:2d}: N={n:8d}  E-block~{n//3:7d}  N^2*8={n*n*8/1e9:7.1f} GB")
