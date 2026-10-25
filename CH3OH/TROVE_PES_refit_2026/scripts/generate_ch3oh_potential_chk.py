#!/usr/bin/env python3
"""generate_ch3oh_potential_chk.py

Exact conversion of a C3v(M)-symmetrised CH3OH PEF (CH3OH_fitpes /
CH3OH_v10 form) into a TROVE potential.chk in the "with_modes" format for
TRANSFORM R-ALPHA132-THETA-S-TAU (modes 1-5 'morse(x-xe)', 6-9 'x-xe',
10-11 'linear', 12 'fourier').

Each symmetrised term (1/6) sum_g prod_m chi_g(m)^p_m T_k(chi_g(12)) is
expanded analytically into monomials of the TROVE coordinates times
cos(k tau)/sin(k tau): methyl permutations act on the exponents, the S_a/S_b
rotation is expanded binomially and the tau shift/inversion via the
cos/sin addition formulae.  No Taylor re-expansion is involved, so the
TROVE potential is identical to the fitted PEF PROVIDED the TROVE input uses
EQUILIBRIUM/SPECPARAM equal to the fitted (rCOe, rOHe, rCHe, aCOHe, aOCHe)
and Morse (a1, a2, b1) parameters (written to a companion file).
"""
import argparse
import math
from collections import defaultdict
import numpy as np
import ch3oh_pes as cp

A_ = 0.5
B_ = 0.5 * math.sqrt(3.0)
P3 = 2 * math.pi / 3
NAMES = ("rCOe", "rOHe", "rCH", "aHOC", "aHCO", "a1", "a2", "b1")

# C3v(M)-2 operations: (r-permutation, S matrix, tau sign, tau shift)
# chi_g(m) = xi(perm[m]) for m in 3..5 and 7..9 (1-based); S' = M S; tau' = s tau + d
OPS = [
    ((3, 4, 5), ((1, 0), (0, 1)), 1, 0.0),
    ((5, 3, 4), ((-A_, -B_), (B_, -A_)), 1, P3),
    ((4, 5, 3), ((-A_, B_), (-B_, -A_)), 1, -P3),
    ((3, 5, 4), ((1, 0), (0, -1)), -1, 0.0),
    ((5, 4, 3), ((-A_, -B_), (-B_, A_)), -1, -P3),
    ((4, 3, 5), ((-A_, B_), (B_, A_)), -1, P3),
]


def binom_expand(c1, c2, n):
  """(c1 Sa + c2 Sb)^n -> {(i, j): coeff} with i+j=n."""
  out = {}
  for j in range(n + 1):
    v = math.comb(n, j) * (c1 ** (n - j)) * (c2 ** j)
    if abs(v) > 1e-300:
      out[(n - j, j)] = v
  return out


def expand_term(p, coeff, out):
  """Add coeff * (1/6) sum_g g.(monomial p) to dict out[(p1..p11,k)]."""
  k = p[11]
  for perm, M, s, d in OPS:
    # rigid permuted exponents: chi_g(3)=xi(perm[0]) -> exponent p[2] goes to perm[0]
    q = [0] * 11
    q[0], q[1], q[5] = p[0], p[1], p[5]
    for slot, src in enumerate(perm):
      q[src - 1] += p[2 + slot]
      q[src + 3] += p[6 + slot]   # bends: modes 7..9 <-> stretches 3..5 (+4)
    # S part
    sa = binom_expand(M[0][0], M[0][1], p[9])
    sb = binom_expand(M[1][0], M[1][1], p[10])
    # tau part: list of (signed k, factor)
    if k >= 0:
      tf = [(k, math.cos(k * d))]
      if k > 0:
        tf.append((-k, -s * math.sin(k * d)))
    else:
      kk = -k
      tf = [(-kk, s * math.cos(kk * d)), (kk, math.sin(kk * d))]
    for (i1, j1), v1 in sa.items():
      for (i2, j2), v2 in sb.items():
        q[9], q[10] = i1 + i2, j1 + j2
        for kt, vt in tf:
          if abs(vt) < 1e-15:
            continue
          key = tuple(q) + (kt,)
          out[key] += coeff * v1 * v2 * vt / 6.0


def read_trove_poten(fname):
  nl, terms = {}, []
  with open(fname) as f:
    inpot = False
    for line in f:
      s = line.split()
      if not s:
        continue
      if s[0].upper() == "POTEN":
        inpot = True
        continue
      if not inpot:
        continue
      if s[0].lower() == "end":
        break
      if s[0] in NAMES:
        nl[s[0]] = float(s[-1])
      elif s[0] in ("f", "h") and len(s) >= 15:
        terms.append((tuple(int(x) for x in s[1:13]), float(s[14])))
  return nl, terms


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("params", help=".trove (POTEN block) or .npz (P, c, nl)")
  ap.add_argument("--out", default="potential.chk")
  ap.add_argument("--thresh", type=float, default=1e-8)
  ap.add_argument("--ncheck", type=int, default=3000)
  a = ap.parse_args()

  if a.params.endswith(".npz"):
    d = np.load(a.params)
    nl = dict(zip(NAMES, d["nl"]))
    terms = [(tuple(int(x) for x in p), float(c)) for p, c in zip(d["P"], d["c"])]
  else:
    nl, terms = read_trove_poten(a.params)
  zero = tuple([0] * 12)
  terms = [(p, c) for p, c in terms if p != zero and c != 0.0]
  print(f"{len(terms)} symmetrised terms; nonlinear: {nl}")

  out = defaultdict(float)
  for p, c in terms:
    expand_term(p, c, out)
  items = [(k, v) for k, v in out.items() if abs(v) >= a.thresh]
  items.sort(key=lambda kv: (sum(kv[0][:11]), kv[0][:11], abs(kv[0][11]), kv[0][11] < 0))
  print(f"expanded into {len(out)} monomials; {len(items)} kept (|c| >= {a.thresh:g})")

  # verification against direct evaluation of the symmetrised PEF
  rng = np.random.default_rng(11)
  xi = np.empty((a.ncheck, 12))
  xi[:, 0:5] = rng.uniform(-0.4, 0.5, (a.ncheck, 5))
  xi[:, 5:9] = rng.uniform(-0.6, 0.6, (a.ncheck, 4))
  xi[:, 9:11] = rng.uniform(-1.0, 1.0, (a.ncheck, 2))
  xi[:, 11] = rng.uniform(0, 2 * np.pi, a.ncheck)
  Ps = np.array([p for p, _ in terms])
  Cs = np.array([c for _, c in terms])
  Vsym = cp.design_matrix(xi, Ps) @ Cs
  Pm = np.array([k for k, _ in items])
  Cm = np.array([v for _, v in items])
  Vmon = cp.design_matrix(xi, Pm, symmetrise=False) @ Cm
  err = np.abs(Vsym - Vmon)
  print(f"check on {a.ncheck} random points: max|V_sym - V_chk| = {err.max():.3e} cm-1"
        f" (|V| up to {np.abs(Vsym).max():.3e})")

  with open(a.out, "w") as f:
    f.write(f"{0:9d}{max(abs(k[11]) for k,_ in items):6d}{len(items):8d}"
            "          <- Npoints,Norder,Ncoeff (header ignored by TROVE)\n")
    for i, (k, v) in enumerate(items, start=1):
      f.write(f"{i:12d}    0  {v:28.18e} " + " ".join(f"{x:4d}" for x in k) + "\n")
    f.write("  987654321       0     0.0000000   " + " ".join("  0" for _ in range(12))
            + " <- End\n")
    f.write(f"{1e-24:14.8E}      <- sparse threshold used\n")
    f.write("End of potential\n")
  with open(a.out + ".trove_coords", "w") as f:
    f.write("EQUILIBRIUM\n")
    f.write(f"re1          0  {nl['rCOe']:16.10f}\n")
    f.write(f"re2          0  {nl['rOHe']:16.10f}\n")
    for i in (3, 4, 5):
      f.write(f"re{i}          0  {nl['rCH']:16.10f}\n")
    f.write(f"alphae1      0  {nl['aHOC']:16.10f} DEG\n")
    for i in (2, 3, 4):
      f.write(f"alphae{i}      0  {nl['aHCO']:16.10f} DEG\n")
    f.write("taue1        0       0.0000      DEG\ntaue2        0     120.0000      DEG\n"
            "taue3        0     240.0000      DEG\nend\n\nSPECPARAM\n")
    f.write(f"beta1        0  {nl['a1']:16.10f}\nbeta2        0  {nl['a2']:16.10f}\n")
    for i in (3, 4, 5):
      f.write(f"beta{i}        0  {nl['b1']:16.10f}\n")
    f.write("END\n")
  print(f"wrote {a.out} and {a.out}.trove_coords")


if __name__ == "__main__":
  main()
