#!/usr/bin/env python3
"""Compare old (962-term) and new fits on the full grid and on the old subset."""
import numpy as np
import ch3oh_pes as cp
import hole_check as hc

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
NB = "/home/user/ExoMol/CH3OH/nb"
crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
E = Eh * cp.HARTREE
Erel = E - E.min()
old_sub = idx <= 19229
# old fit
nl, _, terms = cp.read_params(f"{NB}/00-2D-refitMAPLE-15-09.res")
P = np.array([t[0] for t in terms]); C = np.array([t[2] for t in terms])
Vold = cp.design_matrix(cp.internal_xi(crd, nl, "chart"), P) @ C
rold = E - Vold
rold -= np.median(rold[Erel < 2000])  # absolute-energy offset of the old constant
new = hc.load("../fit/fitW_r0/CH3OH_PES_fitW_r0.trove")
Vnew = new.V_xi(cp.internal_xi(crd, new.nl, "chart"))
rnew = Erel - Vnew
bands = (2000, 5000, 10000, 15000, 20000, 30000, 1e9)
print(f"{'':28s}" + " ".join(f"{b:>8.0f}" for b in bands))
for lab, r, m in (("old fit, old 14608 subset", rold, old_sub),
                  ("new fit, old 14608 subset", rnew, old_sub),
                  ("old fit, all 49537", rold, np.ones_like(old_sub)),
                  ("new fit, all 49537", rnew, np.ones_like(old_sub))):
  out = []
  for b in bands:
    s = m & (Erel < b)
    out.append(np.sqrt(np.mean(r[s] ** 2)))
  print(f"{lab:28s}" + " ".join(f"{x:8.2f}" for x in out))
