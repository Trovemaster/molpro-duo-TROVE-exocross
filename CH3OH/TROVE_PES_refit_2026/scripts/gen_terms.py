#!/usr/bin/env python3
"""gen_terms.py: systematic C3v(M)-symmetrised expansion terms for CH3OH.

Coordinates (TROVE R-ALPHA132-THETA-S-TAU order, 0-based here):
  0 y_CO, 1 y_OH, 2-4 y_CH(1..3), 5 d_COH, 6-8 d_OCH(1..3), 9 S_a, 10 S_b,
  11 tau (Fourier index k: k>=0 cos(k tau), k<0 sin(|k| tau)).

Candidates = canonical methyl patterns x (CO, OH, COH powers) x all S_a/S_b
monomials of a given degree x Fourier functions.  Linear dependencies of the
C3v(M)-symmetrised functions (1/6 sum_g g.f) are removed numerically with a
block pivoted QR on random geometries, keeping the preferred (low-order,
canonical) representatives.
"""
import itertools
import sys
import numpy as np
import scipy.linalg as sla
import ch3oh_pes as cp

GROUP = {0: "CO", 1: "OH", 2: "CH", 3: "CH", 4: "CH", 5: "COH", 6: "OCH",
         7: "OCH", 8: "OCH", 9: "S", 10: "S"}


def canonical_methyl(pr, pb):
  """True if ((pr1,pb1),(pr2,pb2),(pr3,pb3)) is lexicographically maximal
  among the 6 permutations of the three methyl H atoms."""
  pairs = list(zip(pr, pb))
  best = max(itertools.permutations(pairs))
  return tuple(pairs) == best


def ngroups(p):
  return len({GROUP[m] for m in range(11) if p[m] > 0})


def enumerate_candidates(limits, kmax_of_d, kpure=18, dmax=8):
  """limits[n_groups] = max rigid degree allowed for that many groups."""
  cands = []
  # pure torsion
  for k in range(0, kpure + 1):
    cands.append(tuple([0] * 11 + [k]))
  # rigid parts
  for d in range(1, dmax + 1):
    for p in compositions(d, 11):
      ng = ngroups(p)
      if ng not in limits or d > limits[ng]:
        continue
      if not canonical_methyl(p[2:5], p[6:9]):
        continue
      K = kmax_of_d(d, ng, p)
      for k in range(0, K + 1):
        cands.append(tuple(list(p) + [k]))
        if k > 0:
          cands.append(tuple(list(p) + [-k]))
  return cands


def compositions(d, n):
  """All non-negative integer vectors of length n summing to d."""
  if n == 1:
    yield (d,)
    return
  for i in range(d, -1, -1):
    for rest in compositions(d - i, n - 1):
      yield (i,) + rest


def random_xi(n, rng):
  x = np.empty((n, 12))
  x[:, 0:5] = rng.uniform(-0.45, 0.55, (n, 5))
  x[:, 5:9] = rng.uniform(-0.7, 0.7, (n, 4))
  x[:, 9:11] = rng.uniform(-1.2, 1.2, (n, 2))
  x[:, 11] = rng.uniform(0, 2 * np.pi, n)
  return x


def select_independent(cands, npts=14000, tol=1e-9, block=1500, seed=7,
                       verbose=True):
  rng = np.random.default_rng(seed)
  xi = random_xi(npts, rng)
  P = np.array(cands, dtype=int)
  Q = np.zeros((npts, 0))
  keep = []
  for b0 in range(0, len(P), block):
    Pb = P[b0:b0 + block]
    A = cp.design_matrix(xi, Pb)
    nrm = np.linalg.norm(A, axis=0)
    zero = nrm < 1e-12 * np.sqrt(npts)
    A[:, ~zero] /= nrm[~zero]
    if Q.shape[1]:
      A -= Q @ (Q.T @ A)
      A -= Q @ (Q.T @ A)  # re-orthogonalise
    A[:, zero] = 0.0
    # sequential (preference-ordered) Gram-Schmidt inside the block
    accepted = []
    Qb = np.zeros((npts, 0))
    for j in range(A.shape[1]):
      if zero[j]:
        continue
      v = A[:, j].copy()
      if Qb.shape[1]:
        v -= Qb @ (Qb.T @ v)
        v -= Qb @ (Qb.T @ v)
      nv = np.linalg.norm(v)
      if nv > tol * 1e3:   # relative to unit-normalised column
        Qb = np.hstack([Qb, (v / nv)[:, None]])
        accepted.append(b0 + j)
    Q = np.hstack([Q, Qb])
    keep.extend(accepted)
    if verbose:
      print(f"  block {b0:6d}-{b0+len(Pb):6d}: accepted {len(accepted):5d}"
            f"  total {len(keep):6d}", flush=True)
  return [cands[i] for i in keep]


def write_terms(fname, terms):
  with open(fname, "w") as f:
    for p in terms:
      f.write(" f " + " ".join(f"{x:3d}" for x in p) + "\n")


def read_terms(fname):
  out = []
  for l in open(fname):
    s = l.split()
    if s and s[0] == "f":
      out.append(tuple(int(x) for x in s[1:13]))
  return out


if __name__ == "__main__":
  import argparse
  ap = argparse.ArgumentParser()
  ap.add_argument("--l1", type=int, default=8)
  ap.add_argument("--l2", type=int, default=6)
  ap.add_argument("--l3", type=int, default=4)
  ap.add_argument("--l4", type=int, default=4)
  ap.add_argument("--l5", type=int, default=0)
  ap.add_argument("--kmax", default="1:6,2:6,3:4,4:3,5:3,6:3,7:3,8:3")
  ap.add_argument("--kpure", type=int, default=18)
  ap.add_argument("--count-only", action="store_true")
  ap.add_argument("--npts", type=int, default=14000)
  ap.add_argument("--out", default="terms.txt")
  a = ap.parse_args()
  limits = {1: a.l1, 2: a.l2, 3: a.l3, 4: a.l4}
  if a.l5:
    limits[5] = a.l5
  kmap = {int(x.split(":")[0]): int(x.split(":")[1]) for x in a.kmax.split(",")}
  cands = enumerate_candidates(limits, lambda d, ng, p: kmap.get(d, 0),
                               kpure=a.kpure, dmax=max(limits.values()))
  print(f"candidates: {len(cands)}")
  if a.count_only:
    sys.exit()
  terms = select_independent(cands, npts=a.npts)
  print(f"independent symmetrised terms: {len(terms)}")
  write_terms(a.out, terms)
