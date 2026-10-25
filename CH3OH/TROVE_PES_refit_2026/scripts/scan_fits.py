#!/usr/bin/env python3
"""scan_fits.py: model selection on the superset design matrix (one build)."""
import time
import numpy as np
import ch3oh_pes as cp
from gen_terms import read_terms, GROUP
from explore_fit import EX, NL0, ps_weights, solve, band_rms, BANDS

crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
Erel = (Eh - Eh.min()) * cp.HARTREE
P = np.array(read_terms("terms_super.txt"))
t0 = time.time()
xi = cp.internal_xi(crd, NL0, "chart")
F = cp.design_matrix(xi, P)
print(f"superset design matrix {F.shape} in {time.time()-t0:.1f}s", flush=True)
ng = np.array([len({GROUP[m] for m in range(11) if p[m] > 0}) for p in P])
deg = P[:, :11].sum(axis=1)
kab = np.abs(P[:, 11])
rng = np.random.default_rng(123)
test = rng.random(len(Erel)) < 0.1


def subset(l1=8, l2=6, l3=5, l4=4, kmap=None, kpure=18):
  lim = {0: 0, 1: l1, 2: l2, 3: l3, 4: l4}
  m = np.array([deg[i] <= lim.get(ng[i], -1) for i in range(len(P))])
  if kmap:
    m &= np.array([(deg[i] == 0 and kab[i] <= kpure) or
                   (deg[i] > 0 and kab[i] <= kmap.get(deg[i], 0)) for i in range(len(P))])
  return m


def run(label, mask, gamma, emax, ridge):
  t0 = time.time()
  c, info = solve(F[:, mask], Erel, ps_weights(Erel, gamma, emax) * (~test), ridge=ridge)
  res = Erel - F[:, mask] @ c
  tr = band_rms(res, Erel, ~test)
  te = band_rms(res, Erel, test)
  print(f"{label:44s} n={mask.sum():5d} | " + " ".join(f"{x:7.2f}" for x in tr[:6]) +
        " || " + " ".join(f"{x:7.2f}" for x in te[:6]) + f"  ({time.time()-t0:.0f}s)",
        flush=True)
  return c, res


print(f"{'':44s} {'':7s} | train rms: " + " ".join(f"{b:>7.0f}" for b in BANDS[:6]) +
      " || test")
K1 = {1: 9, 2: 6, 3: 4, 4: 3, 5: 3, 6: 3, 7: 2, 8: 2}
K2 = {1: 6, 2: 4, 3: 3, 4: 3, 5: 2, 6: 2, 7: 1, 8: 1}
configs = [
    ("A l8/6/5/4 K1", dict(l1=8, l2=6, l3=5, l4=4, kmap=K1)),
    ("B l8/6/4/4 K1", dict(l1=8, l2=6, l3=4, l4=4, kmap=K1)),
    ("C l8/5/4/4 K1", dict(l1=8, l2=5, l3=4, l4=4, kmap=K1)),
    ("D l8/6/4/0 K1", dict(l1=8, l2=6, l3=4, l4=0, kmap=K1)),
    ("E l8/6/5/4 K2", dict(l1=8, l2=6, l3=5, l4=4, kmap=K2)),
    ("F l8/5/5/4 K1", dict(l1=8, l2=5, l3=5, l4=4, kmap=K1)),
    ("G l6/6/5/4 K1", dict(l1=6, l2=6, l3=5, l4=4, kmap=K1)),
]
for gamma, emax in ((0.0005, 15000), (0.0003, 20000)):
  for ridge in (1e-7, 1e-5):
    print(f"--- PS gamma={gamma} Emax={emax} ridge={ridge}")
    for lab, cfg in configs:
      run(lab, subset(**cfg), gamma, emax, ridge)
