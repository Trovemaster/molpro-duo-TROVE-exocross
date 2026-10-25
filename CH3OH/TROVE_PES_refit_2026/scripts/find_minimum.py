#!/usr/bin/env python3
"""find_minimum.py: equilibrium (global minimum) of a fitted CH3OH PEF."""
import sys
import numpy as np
from scipy.optimize import minimize
import ch3oh_pes as cp
import hole_check as hc

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
pef = hc.load(sys.argv[1])
crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
x0 = cp.internal_xi(crd[np.argmin(Eh)][None, :], pef.nl, "chart")[0]
h = 1e-5


def f_and_g(x):
  pts = np.vstack([x] + [x + h * e for e in np.eye(12)] + [x - h * e for e in np.eye(12)])
  v = pef.V_xi(pts)
  g = (v[1:13] - v[13:25]) / (2 * h)
  return v[0], g


r = minimize(f_and_g, x0, jac=True, method="BFGS", options={"gtol": 1e-6, "maxiter": 500})
x = r.x
nl = pef.nl
q = [nl["rCOe"] - np.log(1 - x[0]) / nl["a1"], nl["rOHe"] - np.log(1 - x[1]) / nl["a2"]]
q += [nl["rCH"] - np.log(1 - x[m]) / nl["b1"] for m in (2, 3, 4)]
q += [(x[5] + nl["aHOC"] * cp.DEG) / cp.DEG]
q += [(x[m] + nl["aHCO"] * cp.DEG) / cp.DEG for m in (6, 7, 8)]
q += [x[9] / cp.DEG, x[10] / cp.DEG, np.mod(x[11], 2 * np.pi) / cp.DEG]
print(f"minimisation: success={r.success} nit={r.nit} V_min = {r.fun:.4f} cm-1"
      f" (relative to the lowest ab initio point)")
names = ["rCO", "rOH", "rCH1", "rCH2", "rCH3", "aCOH", "aOCH1", "aOCH2", "aOCH3",
         "Sa(deg)", "Sb(deg)", "tau(deg)"]
for n_, v in zip(names, q):
  print(f"  {n_:9s} {v:12.6f}")
