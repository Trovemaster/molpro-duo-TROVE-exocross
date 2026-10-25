#!/usr/bin/env python3
"""Count TROVE contracted product functions for CH3OH class structure vs polyad."""
import itertools, math
from collections import Counter

def class_polys(nmodes, vmax, w):
    """polyad values (weighted, real) of class functions: sum v <= vmax"""
    out = []
    for v in itertools.product(range(vmax+1), repeat=nmodes):
        if sum(v) <= vmax:
            out.append(w*sum(v))
    return out

def count(weights, ranges, Pmax, tors_w=0.1):
    # classes: CO(1), OH(1), CH3(3), COH(1), OCH(3), S(2), tors(1, 3 funcs)
    cls = [class_polys(1, ranges['CO'], weights['CO']),
           class_polys(1, ranges['OH'], weights['OH']),
           class_polys(3, ranges['CH'], weights['CH']),
           class_polys(1, ranges['COH'], weights['COH']),
           class_polys(3, ranges['OCH'], weights['OCH']),
           class_polys(2, ranges['S'], weights['S']),
           [tors_w*v for v in range(3)]]
    # convolution
    dist = Counter({0.0: 1})
    for c in cls:
        nd = Counter()
        for p0, n0 in dist.items():
            for p in c:
                nd[round(p0+p, 6)] += n0
        dist = nd
    return sum(n for p, n in dist.items() if int(p+1e-9) <= Pmax), sum(dist.values())

base_r = dict(CO=2, OH=2, CH=2, COH=2, OCH=2, S=4)
for wname, w in [("all 1", dict(CO=1, OH=1, CH=1, COH=1, OCH=1, S=1)),
                 ("2 for OH/CH", dict(CO=1, OH=2, CH=2, COH=1, OCH=1, S=1))]:
    for Pmax in range(2, 9):
        n, tot = count(w, base_r, Pmax)
        print(f"{wname:12s} Pmax={Pmax:2d}: N={n:8d} (full {tot})  ~E block {n//3}")
# frozen-bend reference
r0 = dict(base_r); r0['COH'] = 0; r0['OCH'] = 0
print("frozen bends, all 1, Pmax 23:", count(dict(CO=1, OH=1, CH=1, COH=1, OCH=1, S=1), r0, 23))
