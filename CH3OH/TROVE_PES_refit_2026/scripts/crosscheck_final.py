#!/usr/bin/env python3
"""Cross-check: final .trove parameters (python + fortran evaluator) vs 00.res calc."""
import gzip
import numpy as np
import ch3oh_pes as cp
import hole_check as hc

T = "../fit/fitW40_r2/CH3OH_PES_fitW40_r2.trove"
R = "../fit/fitW40_r2/00.res"
pef = hc.load(T)
rows = []
with open(R) as f:
  for l in f:
    s = l.split()
    if len(s) == 16:
      try:
        rows.append([float(x) for x in s])
      except ValueError:
        pass
r = np.array(rows)[:49537]
crd = r[:, :12]
calc = r[:, 13]
xi = cp.internal_xi(crd, pef.nl, "chart")
idx = np.random.default_rng(0).choice(len(r), 3000, replace=False)
Vf = pef.V_xi(xi[idx])
Vp = cp.design_matrix(xi[idx], pef.P) @ pef.c
print("fortran vs python evaluator: max diff", np.abs(Vf - Vp).max())
print("evaluator vs 00.res calc (coords rounded to 1e-6/1e-6 deg): max diff",
      np.abs(Vf - calc[idx]).max(), " median", np.median(np.abs(Vf - calc[idx])))
