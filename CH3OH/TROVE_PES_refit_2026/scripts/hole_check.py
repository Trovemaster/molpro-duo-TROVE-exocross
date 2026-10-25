#!/usr/bin/env python3
"""hole_check.py: large-displacement / hole diagnostics of a fitted CH3OH PEF.

Checks (V relative to the lowest ab initio energy, cm-1):
 1. 1D cuts along all 12 coordinates through the staggered reference
    (tau = 60 deg) over the TROVE basis borders.
 2. 2D cuts for all coordinate pairs on moderate ranges.
 3. Random multi-dimensional sampling inside a "physical" region defined by
    the separable (sum of 1D cuts) energy V_sep < Ecap; reports points with
    V < Vlow (default -10 cm-1) or V < frac*V_sep.
 4. Batched projected-gradient local minimisation from random starts inside
    the box, reporting any minimum below Vlow.
"""
import argparse
import os
import subprocess
import tempfile
import numpy as np
import ch3oh_pes as cp

DEG = cp.DEG
EVAL_EXE = "/home/user/ExoMol/CH3OH/pes_refit/fit/pes_eval.x"


class PEF:
  def __init__(self, P, c, nl, trove_file=None):
    self.P = np.asarray(P, dtype=int)
    self.c = np.asarray(c, dtype=float)
    self.nl = nl
    self.tmp = tempfile.mkdtemp(prefix="pes_eval_")
    if trove_file is None:
      trove_file = os.path.join(self.tmp, "pef.trove")
      with open(trove_file, "w") as f:
        f.write("POTEN\n")
        for p, v in zip(self.P, self.c):
          f.write(" f " + " ".join(f"{x:3d}" for x in p) + f"  1  {v:.17e}\n")
        f.write("end\n")
    self.trove_file = trove_file

  def V_xi(self, xi, chunk=None):
    xi = np.ascontiguousarray(xi, dtype=np.float64)
    fx = os.path.join(self.tmp, "xi.bin")
    fv = os.path.join(self.tmp, "v.bin")
    with open(fx, "wb") as f:
      np.array([xi.shape[0]], dtype=np.int64).tofile(f)
      xi.tofile(f)   # row-major (N,12) == Fortran column-major (12,N)
    subprocess.run([EVAL_EXE, self.trove_file, fx, fv], check=True)
    return np.fromfile(fv, dtype=np.float64, count=xi.shape[0])

  # physical <-> fit coordinates
  def phys_to_xi(self, q):
    """q: rCO, rOH, r1, r2, r3 (A), aCOH, b1, b2, b3 (rad), Sa, Sb, tau (rad)
    with r2,r3,b2,b3 in the fit (TROVE) labelling."""
    nl = self.nl
    xi = np.empty_like(q)
    xi[:, 0] = 1 - np.exp(-nl["a1"] * (q[:, 0] - nl["rCOe"]))
    xi[:, 1] = 1 - np.exp(-nl["a2"] * (q[:, 1] - nl["rOHe"]))
    for m in (2, 3, 4):
      xi[:, m] = 1 - np.exp(-nl["b1"] * (q[:, m] - nl["rCH"]))
    xi[:, 5] = q[:, 5] - nl["aHOC"] * DEG
    for m in (6, 7, 8):
      xi[:, m] = q[:, m] - nl["aHCO"] * DEG
    xi[:, 9:12] = q[:, 9:12]
    return xi

  def ref_q(self):
    nl = self.nl
    return np.array([nl["rCOe"], nl["rOHe"], nl["rCH"], nl["rCH"], nl["rCH"],
                     nl["aHOC"] * DEG, nl["aHCO"] * DEG, nl["aHCO"] * DEG,
                     nl["aHCO"] * DEG, 0.0, 0.0, 60 * DEG])


NAMES = ["rCO", "rOH", "rCH1", "rCH2", "rCH3", "aCOH", "aOCH1", "aOCH2",
         "aOCH3", "Sa", "Sb", "tau"]
# TROVE basis borders (CH3OH working input) and S/tau ranges
BORDERS = [(1.0, 3.6), (0.55, 3.10), (0.55, 3.3), (0.55, 3.3), (0.55, 3.3),
           (60 * DEG, 170 * DEG), (60 * DEG, 160 * DEG), (60 * DEG, 160 * DEG),
           (60 * DEG, 160 * DEG), (-100 * DEG, 100 * DEG),
           (-100 * DEG, 100 * DEG), (0.0, 2 * np.pi)]
# box for multi-D sampling / minimisation
BOX = [(1.15, 2.3), (0.70, 1.80), (0.80, 1.90), (0.80, 1.90), (0.80, 1.90),
       (70 * DEG, 165 * DEG), (75 * DEG, 155 * DEG), (75 * DEG, 155 * DEG),
       (75 * DEG, 155 * DEG), (-70 * DEG, 70 * DEG), (-70 * DEG, 70 * DEG),
       (0.0, 2 * np.pi)]


def cuts_1d(pef, npts=241, verbose=True):
  q0 = pef.ref_q()
  V0 = pef.V_xi(pef.phys_to_xi(q0[None, :]))[0]
  report = {}
  for m in range(12):
    lo, hi = BORDERS[m]
    g = np.linspace(lo, hi, npts)
    q = np.tile(q0, (npts, 1))
    q[:, m] = g
    V = pef.V_xi(pef.phys_to_xi(q))
    # local minima (excluding ends)
    lm = np.where((V[1:-1] < V[:-2]) & (V[1:-1] < V[2:]))[0] + 1
    ends = (V[0], V[-1])
    report[NAMES[m]] = (V.min(), g[np.argmin(V)], len(lm), ends)
    if verbose:
      unit = 1.0 if m < 5 else 1 / DEG
      s = ", ".join(f"{g[i]*unit:.3f}:{V[i]:.0f}" for i in lm[:4])
      print(f"  {NAMES[m]:6s} Vmin={V.min():10.1f} at {g[np.argmin(V)]*unit:8.3f}"
            f"  ends=({ends[0]:10.0f},{ends[1]:10.0f})  local minima: {s}")
  return V0, report


def cuts_2d(pef, n=41, frac=0.8, verbose=True):
  q0 = pef.ref_q()
  worst = []
  for i in range(12):
    for j in range(i + 1, 12):
      gi = np.linspace(*shrink(BOX[i], frac), n)
      gj = np.linspace(*shrink(BOX[j], frac), n)
      Gi, Gj = np.meshgrid(gi, gj, indexing="ij")
      q = np.tile(q0, (n * n, 1))
      q[:, i] = Gi.ravel()
      q[:, j] = Gj.ravel()
      V = pef.V_xi(pef.phys_to_xi(q))
      worst.append((V.min(), NAMES[i], NAMES[j]))
  worst.sort()
  if verbose:
    print("  lowest 2D-cut values (pairs):")
    for v, a, b in worst[:8]:
      print(f"    {a:6s}-{b:6s}: Vmin = {v:10.1f}")
  return worst


def shrink(b, f):
  c = 0.5 * (b[0] + b[1])
  h = 0.5 * (b[1] - b[0]) * f
  return (c - h, c + h)


def sample_box(n, rng, box=BOX):
  q = np.empty((n, 12))
  for m in range(12):
    q[:, m] = rng.uniform(box[m][0], box[m][1], n)
  return q


def separable_energy(pef, q, q0, V0):
  """Sum of 1D-cut energies (relative to V0)."""
  Vs = np.zeros(q.shape[0])
  for m in range(11):
    qq = np.tile(q0, (q.shape[0], 1))
    qq[:, m] = q[:, m]
    qq[:, 11] = q[:, 11]
    Vs += pef.V_xi(pef.phys_to_xi(qq)) - V0
  return Vs


def random_check(pef, n=200000, ecap=30000.0, vlow=-10.0, seed=3, verbose=True):
  rng = np.random.default_rng(seed)
  q0 = pef.ref_q()
  V0 = pef.V_xi(pef.phys_to_xi(q0[None, :]))[0]
  # sample with Gaussian displacements (scaled to the box) around q0
  q = np.tile(q0, (n, 1))
  for m in range(11):
    lo, hi = BOX[m]
    sig = 0.35 * (hi - lo)
    q[:, m] = np.clip(q0[m] + rng.normal(0, sig, n) * rng.random(n), lo, hi)
  q[:, 11] = rng.uniform(0, 2 * np.pi, n)
  V = pef.V_xi(pef.phys_to_xi(q))
  Vs = separable_energy(pef, q, q0, V0) + V0
  phys = Vs < ecap
  nlow = np.sum(phys & (V < vlow))
  nrel = np.sum(phys & (V < 0.3 * Vs) & (Vs > 5000))
  if verbose:
    print(f"  random: {phys.sum()} of {n} points with V_sep<{ecap:.0f};"
          f" V<{vlow}: {nlow};  V<0.3*V_sep (V_sep>5000): {nrel};"
          f"  min V = {V[phys].min():.1f}")
  return V[phys].min(), nlow, nrel


def minimise_from_random(pef, n=2000, steps=150, seed=5, vlow=-10.0, verbose=True):
  """Projected gradient descent (finite differences) in the box."""
  rng = np.random.default_rng(seed)
  q = sample_box(n, rng)
  lo = np.array([b[0] for b in BOX])
  hi = np.array([b[1] for b in BOX])
  h = np.array([1e-3] * 5 + [1e-3] * 6 + [1e-3])
  step = np.array([0.02] * 5 + [0.02] * 6 + [0.05])
  V = pef.V_xi(pef.phys_to_xi(q))
  for it in range(steps):
    G = np.zeros_like(q)
    for m in range(12):
      qp = q.copy(); qp[:, m] += h[m]
      qm = q.copy(); qm[:, m] -= h[m]
      G[:, m] = (pef.V_xi(pef.phys_to_xi(qp)) - pef.V_xi(pef.phys_to_xi(qm))) / (2 * h[m])
    gn = np.sqrt(np.sum((G * step) ** 2, axis=1)) + 1e-30
    qn = q - (G * step) / gn[:, None] * step.max() * 0.5
    qn = np.clip(qn, lo, hi)
    Vn = pef.V_xi(pef.phys_to_xi(qn))
    better = Vn < V
    q[better] = qn[better]
    V[better] = Vn[better]
    step = step * np.where(better.mean() > 0.5, 1.05, 0.8)
  bad = V < vlow
  if verbose:
    print(f"  minimisation: {n} starts, lowest V found = {V.min():.1f};"
          f" minima below {vlow}: {bad.sum()}")
    if bad.any():
      i = np.argmin(V)
      qq = q[i].copy()
      qq[5:] = qq[5:] / DEG
      print("   worst point:", np.array2string(qq, precision=3))
  return V.min(), bad.sum()


def load(fname):
  keys = ["rCOe", "rOHe", "rCH", "aHOC", "aHCO", "a1", "a2", "b1"]
  if fname.endswith(".npz"):
    d = np.load(fname)
    return PEF(d["P"], d["c"], dict(zip(keys, d["nl"])))
  # .trove file written by CH3OH_fitpes
  nl, P, C, v0 = {}, [], [], 0.0
  inpot = False
  for line in open(fname):
    s = line.split()
    if not s:
      continue
    if line.startswith("! Constant term V0"):
      v0 = float(line.split("=")[1].split()[0])
    if s[0].upper() == "POTEN":
      inpot = True
      continue
    if not inpot:
      continue
    if s[0].lower() == "end":
      break
    if s[0] in keys:
      nl[s[0]] = float(s[-1])
    elif s[0] in ("f", "h"):
      p = [int(x) for x in s[1:13]]
      v = float(s[14])
      if all(x == 0 for x in p) and v == 0.0:
        v = v0
      P.append(p)
      C.append(v)
  return PEF(np.array(P), np.array(C), nl)


if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument("npz")
  ap.add_argument("--no-min", action="store_true")
  ap.add_argument("--nmin", type=int, default=600)
  ap.add_argument("--steps", type=int, default=80)
  ap.add_argument("--nrand", type=int, default=40000)
  a = ap.parse_args()
  pef = load(a.npz)
  print("1D cuts through the staggered reference (TROVE borders):")
  cuts_1d(pef)
  cuts_2d(pef)
  random_check(pef, n=a.nrand)
  if not a.no_min:
    minimise_from_random(pef, n=a.nmin, steps=a.steps)
