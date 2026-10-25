#!/usr/bin/env python3
"""make_pseudo.py: pseudo-points that keep the PEF hole-free where no ab
initio data exist.

Pseudo-points are sampled inside the physical shell V_sep < Ecap (V_sep =
separable energy from the 1D cuts of the current fit, which are well
determined by the ab initio data), restricted to the torsional fundamental
domain tau in [60,120] deg used by the grid, and kept only if they are far
from every ab initio point (normalised coordinates).  Hole geometries found
by shell_holes.py are added (mapped into the fundamental domain by the
C3v(M) operations).  Target energy: E = ratio * V_sep (ratio = 0.9 by
default; on the grid E_ab/V_sep has median 1.0 and minimum 0.52).

Output: data lines in the grid (Molpro z-matrix) format, energy in cm-1
relative to E_min, with the absolute weight in column 14.
"""
import argparse
import numpy as np
from scipy.spatial import cKDTree
import ch3oh_pes as cp
import hole_check as hc
import shell_holes as sh

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
DEG = cp.DEG
# normalisation of fit coordinates for distances
SCALE = np.array([0.08, 0.08, 0.08, 0.08, 0.08, 0.10, 0.10, 0.10, 0.10, 0.12, 0.12, 0.25])


def xi_to_grid(pef, xi):
  """Fit coordinates -> Molpro grid columns (file labelling, deg)."""
  nl = pef.nl
  q = np.empty((xi.shape[0], 12))
  rCO = nl["rCOe"] - np.log(1 - xi[:, 0]) / nl["a1"]
  rOH = nl["rOHe"] - np.log(1 - xi[:, 1]) / nl["a2"]
  r = [nl["rCH"] - np.log(1 - xi[:, m]) / nl["b1"] for m in (2, 3, 4)]
  b = [(xi[:, m] + nl["aHCO"] * DEG) / DEG for m in (6, 7, 8)]
  t1, t2, t3 = cp.chart_to_dihedrals(xi[:, 9], xi[:, 10], xi[:, 11])
  q[:, 0], q[:, 1] = rCO, rOH
  q[:, 2], q[:, 3], q[:, 4] = r[0], r[2], r[1]          # file r2 = code r3
  q[:, 5] = (xi[:, 5] + nl["aHOC"] * DEG) / DEG
  q[:, 6], q[:, 7], q[:, 8] = b[0], b[2], b[1]          # file b2 = code b3
  q[:, 9], q[:, 10], q[:, 11] = t1 / DEG, t3 / DEG, t2 / DEG  # file d2=code th3
  return q


def to_fundamental(xi):
  """Map tau into [60,120] deg with the C3v(M) operations."""
  ops = cp.sym_ops_xi(xi)
  out = xi.copy()
  done = np.zeros(len(xi), bool)
  for o in ops:
    t = np.mod(o[:, 11], 2 * np.pi)
    ok = (~done) & (t >= 60 * DEG - 1e-9) & (t <= 120 * DEG + 1e-9)
    out[ok] = o[ok]
    out[ok, 11] = t[ok]
    done |= ok
  return out, done


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("pef")
  ap.add_argument("--n", type=int, default=600000, help="random samples")
  ap.add_argument("--ecap", type=float, default=30000.0)
  ap.add_argument("--emin-sep", type=float, default=4000.0)
  ap.add_argument("--dmin", type=float, default=1.5, help="min normalised distance to data")
  ap.add_argument("--ratio", type=float, default=0.9)
  ap.add_argument("--weight", type=float, default=0.02)
  ap.add_argument("--holes", nargs="*", default=[])
  ap.add_argument("--hole-weight", type=float, default=0.005)
  ap.add_argument("--max-keep", type=int, default=0)
  ap.add_argument("--seed", type=int, default=11)
  ap.add_argument("--out", required=True)
  a = ap.parse_args()
  pef = hc.load(a.pef)
  rng = np.random.default_rng(a.seed)
  q0 = pef.ref_q()
  V0 = pef.V_xi(pef.phys_to_xi(q0[None, :]))[0]

  crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
  xd = cp.internal_xi(crd, pef.nl, "chart")
  tree = cKDTree(xd / SCALE)

  q = sh.sample(pef, a.n, rng, q0)
  q[:, 11] = rng.uniform(60 * DEG, 120 * DEG, a.n)
  Vs = sh.vsep(pef, q, q0, V0)
  sel = (Vs < a.ecap) & (Vs > a.emin_sep)
  xs = pef.phys_to_xi(q[sel])
  d, _ = tree.query(xs / SCALE, k=1)
  far = d > a.dmin
  xs, Vss = xs[far], Vs[sel][far]
  if len(xs) > a.max_keep:
    a.max_keep = max(a.max_keep, 0)
    k = rng.choice(len(xs), a.max_keep, replace=False)
    xs, Vss = xs[k], Vss[k]
  w = np.full(len(xs), a.weight)
  print(f"random pseudo-points: {sel.sum()} in shell, {far.sum()} far from data, kept {len(xs)}")
  allx, allE, allw = [xs], [a.ratio * Vss], [w]
  for hf in a.holes:
    h = np.load(hf)
    if len(h["q"]) == 0:
      continue
    xh = pef.phys_to_xi(h["q"])
    xh, ok = to_fundamental(xh)
    xh = xh[ok]
    Vh = h["Vsep"][ok]
    m = (Vh > 1500) & (Vh < 2 * a.ecap)
    allx.append(xh[m]); allE.append(a.ratio * np.minimum(Vh[m], a.ecap * 1.5))
    allw.append(np.full(m.sum(), a.hole_weight))
    print(f"hole points from {hf}: {m.sum()}")
  X = np.concatenate(allx); E = np.concatenate(allE); W = np.concatenate(allw)
  G = xi_to_grid(pef, X)
  # consistency check: back-transform
  xb = cp.internal_xi(G, pef.nl, "chart")
  dd = np.abs(np.mod(xb - X + np.pi, 2 * np.pi) - np.pi)
  print(f"round-trip check max |dxi| = {dd.max():.2e}")
  with open(a.out, "w") as f:
    for g, e, ww in zip(G, E, W):
      f.write(" ".join(f"{x:.8f}" for x in g) + f" {e:16.6f} {ww:.6f}  -1\n")
  print(f"wrote {len(G)} pseudo-points to {a.out} (E range {E.min():.0f}-{E.max():.0f})")


if __name__ == "__main__":
  main()
