#!/usr/bin/env python3
"""Step-1 consistency tests for the CH3OH PES coordinates.

T1  Reproduce the 00-2D-refitMAPLE-15-09.res calc values with the python
    re-implementation of CH3OH_v10 (diff_V_tau + XY3_IV).
T2  Invariance of the fitted PES under the actual C3v(M) permutation-
    inversion operations applied to the ab initio grid geometries
    (raw-tau mapping of CH3OH_v10 vs chart-consistent tau).
T3  Census of grid points whose raw tau differs from the TROVE chart tau.
T4  Working potential.chk vs the fit: (a) with the fit's own Morse/eq.
    parameters, (b) with the TROVE EQUILIBRIUM/SPECPARAM of the working
    TROVE input.
"""
import sys
import numpy as np
import ch3oh_pes as cp

NB = "/home/user/ExoMol/CH3OH/nb"
EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"

nl, _, terms = cp.read_params(f"{NB}/00-2D-refitMAPLE-15-09.res")
print("nonlinear params (.res):", nl)
P = np.array([t[0] for t in terms])
C = np.array([t[2] for t in terms])
ve = C[np.all(P == 0, axis=1)][0]
mask = ~np.all(P == 0, axis=1)
P, C = P[mask], C[mask]
print(f"{len(C)} linear terms (+ constant ve={ve})")

# ---- the 14608 points used in that fit (same order as the .res) --------
inp = f"{NB}/CH3OH-2D-NewKin-TorsionSwap.inp"
rows = []
with open(inp) as f:
  lines = f.readlines()
for l in lines[1063:]:
  s = l.split()
  if len(s) >= 14 and float(s[0]) > 0:
    rows.append([float(x) for x in s[:14]])
rows = np.array(rows)
crd, Eh = rows[:, :12], rows[:, 12]
res = []
with open(f"{NB}/00-2D-refitMAPLE-15-09.res") as f:
  for l in f:
    s = l.split()
    if len(s) == 16:
      try:
        res.append([float(x) for x in s])
      except ValueError:
        pass
res = np.array(res)
print("points: inp", len(crd), " res", len(res))
n = min(len(crd), len(res))
crd, Eh, res = crd[:n], Eh[:n], res[:n]
assert np.allclose(crd[:, :5], res[:, :5], atol=2e-6)

# ---- T1 -----------------------------------------------------------------
xi = cp.internal_xi(crd, nl, "raw")
F = cp.design_matrix(xi, P)
V = F @ C
calc = res[:, 13]
obs = res[:, 12]
erel = obs - obs.min()
d = V - calc
sel = erel < 20000
print("\nT1: |V_python - calc_res| max/rms (all)       = %.3e / %.3e cm-1"
      % (np.abs(d).max(), np.sqrt(np.mean(d**2))))
print("T1: |V_python - calc_res| max/rms (E<20000)   = %.3e / %.3e cm-1"
      % (np.abs(d[sel]).max(), np.sqrt(np.mean(d[sel]**2))))


# ---- T2/T3 on the full 49537-point grid -----------------------------------
crdA, EhA, idxA = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
xr = cp.internal_xi(crdA, nl, "raw")
xc = cp.internal_xi(crdA, nl, "chart")
dt = np.mod(xr[:, 11] - xc[:, 11] + np.pi, 2 * np.pi) - np.pi
nbad = np.sum(np.abs(dt) > 1e-6)
print(f"\nT3: grid points with raw tau != chart tau: {nbad} of {len(dt)}")
if nbad:
  vals, cnts = np.unique(np.round(dt[np.abs(dt) > 1e-6] / cp.DEG), return_counts=True)
  print("    offsets (deg):", dict(zip(vals, cnts)))
  Erel = (EhA - EhA.min()) * cp.HARTREE
  bad = np.abs(dt) > 1e-6
  print("    E_rel of affected points: min %.1f  median %.1f cm-1"
        % (Erel[bad].min(), np.median(Erel[bad])))
  th1 = crdA[bad, 9]
  print("    theta1 range of affected points: %.2f .. %.2f deg"
        % (th1.min(), th1.max()))


def apply_op(c, op):
  """Actual C3v(M) operation on Molpro grid coordinates (file labels)."""
  c = c.copy()
  r = c[:, 2:5].copy(); b = c[:, 6:9].copy(); t = c[:, 9:12].copy()
  perms = {"E": (0, 1, 2), "(123)": (2, 0, 1), "(132)": (1, 2, 0),
           "(12)*": (1, 0, 2), "(13)*": (2, 1, 0), "(23)*": (0, 2, 1)}
  p = perms[op]
  c[:, 2:5] = r[:, p]; c[:, 6:9] = b[:, p]; c[:, 9:12] = t[:, p]
  if op.endswith("*"):
    c[:, 9:12] = np.mod(-c[:, 9:12], 360.0)
  return c


rng = np.random.default_rng(1)
sub = rng.choice(len(crdA), 3000, replace=False)
c0 = crdA[sub]
Erel0 = (EhA[sub] - EhA.min()) * cp.HARTREE
for mode in ("raw", "chart"):
  V0 = cp.design_matrix(cp.internal_xi(c0, nl, mode), P) @ C
  print(f"\nT2 [{mode} tau]: max |V(g.x)-V(x)| over 3000 random grid points")
  for op in ("(123)", "(132)", "(12)*", "(13)*", "(23)*"):
    Vg = cp.design_matrix(cp.internal_xi(apply_op(c0, op), nl, mode), P) @ C
    dv = np.abs(Vg - V0)
    lo = Erel0 < 15000
    print(f"   {op:6s}: all {dv.max():12.4e}   E<15000: {dv[lo].max():12.4e}"
          f"   n(>1e-4) = {np.sum(dv > 1e-4)}")

# ---- T4 -------------------------------------------------------------------
Pc, Cc = cp.read_potential_chk(f"{NB}/potential.chk")
print(f"\nT4: potential.chk has {len(Cc)} terms; max rigid power "
      f"{Pc[:, :11].max()}, fourier k range {Pc[:, 11].min()}..{Pc[:, 11].max()}")
eq_trove = {"re1": 1.42194259, "re2": 0.95640043, "re3": 1.09066795,
            "alphae1": 108.09748505, "alphae2": 110.28648298}
beta_trove = (2.16893551, 2.57276067, 1.78266276)
eq_fit = {"re1": nl["rCOe"], "re2": nl["rOHe"], "re3": nl["rCH"],
          "alphae1": nl["aHOC"], "alphae2": nl["aHCO"]}
beta_fit = (nl["a1"], nl["a2"], nl["b1"])
for lab, eq, bt in (("fit eq/Morse", eq_fit, beta_fit),
                    ("TROVE input EQUILIBRIUM/SPECPARAM", eq_trove, beta_trove)):
  chi = cp.trove_chi_pot(crd, eq, bt, tau_mode="chart")
  Vchk = cp.eval_potential_chk(chi, Pc, Cc)
  dd = Vchk - V
  for emax in (5000, 10000, 20000):
    s = erel < emax
    print(f"   [{lab}] E<{emax:5d}: max|Vchk-Vfit| = {np.abs(dd[s]).max():10.3f}"
          f"  rms = {np.sqrt(np.mean(dd[s]**2)):9.3f} cm-1  (n={s.sum()})")
