#!/usr/bin/env python3
"""Calibrate E_ab / V_sep on the ab initio grid (multi-D points)."""
import sys
import numpy as np
import ch3oh_pes as cp
import hole_check as hc
import shell_holes as sh

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
pef = hc.load(sys.argv[1])
crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
Erel = (Eh - Eh.min()) * cp.HARTREE
xi = cp.internal_xi(crd, pef.nl, "chart")
# physical coordinates q (fit labelling) from xi
nl = pef.nl
q = np.empty_like(xi)
q[:, 0] = nl["rCOe"] - np.log(1 - xi[:, 0]) / nl["a1"]
q[:, 1] = nl["rOHe"] - np.log(1 - xi[:, 1]) / nl["a2"]
for m in (2, 3, 4):
  q[:, m] = nl["rCH"] - np.log(1 - xi[:, m]) / nl["b1"]
q[:, 5] = xi[:, 5] + nl["aHOC"] * cp.DEG
for m in (6, 7, 8):
  q[:, m] = xi[:, m] + nl["aHCO"] * cp.DEG
q[:, 9:12] = xi[:, 9:12]
q0 = pef.ref_q()
V0 = pef.V_xi(pef.phys_to_xi(q0[None, :]))[0]
Vs = sh.vsep(pef, q, q0, V0)
V = pef.V_xi(xi)
m = Vs > 5000
r = Erel[m] / Vs[m]
print(f"multi-D points with V_sep>5000: {m.sum()}")
for p in (1, 5, 10, 25, 50, 75, 90):
  print(f"  percentile {p:3d}: E_ab/V_sep = {np.percentile(r, p):.3f}")
print("  min ratio:", r.min(), " fraction with ratio<0.5:", np.mean(r < 0.5))
rf = V[m] / Vs[m]
print("  fit: min V/V_sep on data =", rf.min())
