#!/usr/bin/env python3
"""count_basis3.py: CH3OH product-basis size with polyad AND zero-order energy cut (approximate class energies)."""
import itertools
from collections import defaultdict

W = dict(CO=1, OH=2, CH=2, COH=1, OCH=1, S=1)
OM = dict(CO=1100., OH=3600., CH=2900., COH=1270., OCH=1400., S=1500.)
NM = dict(CO=1, OH=1, CH=3, COH=1, OCH=3, S=2)

def cls(name, vmax):
  out = []
  for v in itertools.product(range(vmax + 1), repeat=NM[name]):
    s = sum(v)
    if s <= vmax:
      out.append((int(W[name] * s + 1e-9), OM[name] * s))
  return out

def tors(nf, wt, etrip=(0., 250., 450., 700., 1000.)):
  return [(int(wt * n + 1e-9), etrip[n // 3]) for n in range(nf)]

def count(r, nt, wt, pmax, ecut):
  classes = [cls(k, r[k]) for k in ("CO", "OH", "CH", "COH", "OCH", "S")] + [tors(nt, wt)]
  # dynamic programming over (polyad, energy-bin)
  dist = defaultdict(int); dist[(0, 0.0)] = 1
  for c in classes:
    nd = defaultdict(int)
    for (p0, e0), n0 in dist.items():
      for p, e in c:
        if p0 + p <= pmax and e0 + e <= ecut:
          nd[(p0 + p, round(e0 + e, 1))] += n0
    dist = nd
  return sum(dist.values())

r4 = dict(CO=4, OH=4, CH=4, COH=4, OCH=4, S=4)
for nt, wt in ((9, 0.34), (9, 0.1), (12, 0.34)):
  for pmax in (8, 9, 10):
    row = []
    for ecut in (8000, 10000, 12000, 14000, 99999):
      row.append(count(r4, nt, wt, pmax, ecut))
    print(f"tors {nt} w={wt} Pmax={pmax}: N(Ecut=8k,10k,12k,14k,inf) = {row}")
