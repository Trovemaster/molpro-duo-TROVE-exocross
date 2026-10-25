#!/usr/bin/env python3
"""plot_cuts.py: 1D cuts of a fitted PEF vs the ab initio 1D grid points,
extended over the TROVE basis borders (large-displacement behaviour)."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ch3oh_pes as cp
import hole_check as hc

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
pef = hc.load(sys.argv[1])
out = sys.argv[2] if len(sys.argv) > 2 else "cuts_1d.png"
crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
Erel = (Eh - Eh.min()) * cp.HARTREE
xi = cp.internal_xi(crd, pef.nl, "chart")
xr = xi[0]
d = xi - xr
d[:, 11] = np.mod(d[:, 11] + np.pi, 2 * np.pi) - np.pi
labels = ["r(CO) / A", "r(OH) / A", "r(CH1) / A", "r(CH2) / A", "r(CH3) / A",
          "a(COH) / deg", "a(OCH1) / deg", "a(OCH2) / deg", "a(OCH3) / deg",
          "S_a / deg", "S_b / deg", "tau / deg"]
# grid reference in physical coordinates (fit labelling)
q_ref = np.array([crd[0, 0], crd[0, 1], crd[0, 2], crd[0, 4], crd[0, 3],
                  crd[0, 5] * cp.DEG, crd[0, 6] * cp.DEG, crd[0, 8] * cp.DEG,
                  crd[0, 7] * cp.DEG, xr[9], xr[10], xr[11]])
modes = [0, 1, 2, 5, 6, 9, 10, 11]
fig, axs = plt.subplots(2, 4, figsize=(17, 8))
for ax, m in zip(axs.ravel(), modes):
  others = [j for j in range(12) if j != m]
  sel = np.all(np.abs(d[:, others]) < 1e-7, axis=1) & (np.abs(d[:, m]) > 1e-9)
  lo, hi = hc.BORDERS[m]
  g = np.linspace(lo, hi, 600)
  q = np.tile(q_ref, (len(g), 1))
  q[:, m] = g
  V = pef.V_xi(pef.phys_to_xi(q))
  unit = 1.0 if m < 5 else 1.0 / cp.DEG
  ax.plot(g * unit, V, "b-", lw=1.5, label="fit")
  if sel.any():
    if m < 5:
      xq = {0: crd[sel, 0], 1: crd[sel, 1], 2: crd[sel, 2]}[m]
    elif m in (5, 6):
      xq = {5: crd[sel, 5], 6: crd[sel, 6]}[m]
    else:
      xq = (xr[m] + d[sel, m]) / cp.DEG
    ax.plot(xq, Erel[sel], "ro", ms=3.5, label=f"ab initio ({sel.sum()})")
  ax.set_xlabel(labels[m])
  ax.set_ylabel("E / cm$^{-1}$")
  top = max(45000, np.nanmax(Erel[sel]) * 1.15 if sel.any() else 45000)
  ax.set_ylim(-500, min(top, 1.2e5) if m != 11 else 600)
  ax.axhline(0, color="k", lw=0.5)
  ax.legend(fontsize=8)
fig.suptitle("CH3OH PEF: 1D cuts through the grid reference (other coordinates fixed), "
             "TROVE basis borders")
fig.tight_layout()
fig.savefig(out, dpi=110)
print("saved", out)
