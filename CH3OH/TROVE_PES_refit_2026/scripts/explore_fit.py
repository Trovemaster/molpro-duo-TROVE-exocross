#!/usr/bin/env python3
"""explore_fit.py: fast exploratory linear fits (nonlinear params fixed)."""
import argparse
import time
import numpy as np
import scipy.linalg as sla
import ch3oh_pes as cp
from gen_terms import read_terms

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
NL0 = {"rCOe": 1.42167424, "rOHe": 0.95631449, "rCH": 1.09061196,
       "aHOC": 108.10245297, "aHCO": 110.28972802, "a1": 2.17023887,
       "a2": 2.55392176, "b1": 1.85744195}
BANDS = (2000, 5000, 10000, 15000, 20000, 30000, 1e9)


def ps_weights(Erel, gamma, emax):
  if emax <= 0:
    return np.ones_like(Erel)
  return (np.tanh(-gamma * (Erel - emax)) + 1.000200020002) / 2.000200020002


def band_rms(res, Erel, mask=None):
  out = []
  for b in BANDS:
    s = Erel < b
    if mask is not None:
      s &= mask
    out.append(np.sqrt(np.mean(res[s] ** 2)) if s.any() else np.nan)
  return out


def solve(F, y, w, ridge=0.0, cutoff=0.0):
  sw = np.sqrt(w)
  A = F * sw[:, None]
  b = y * sw
  N = A.T @ A
  r = A.T @ b
  d = np.sqrt(np.maximum(np.diag(N), 1e-300))
  Ns = N / np.outer(d, d)
  rs = r / d
  if cutoff > 0:
    lam, V = np.linalg.eigh(Ns)
    keep = lam > cutoff * lam.max()
    c = V[:, keep] @ ((V[:, keep].T @ rs) / lam[keep])
    info = f"eig kept {keep.sum()}/{len(lam)}, cond={lam.max()/lam[lam>0].min():.2e}"
  else:
    Ns[np.diag_indices_from(Ns)] += ridge
    c = sla.cho_solve(sla.cho_factor(Ns, lower=True), rs)
    info = "cholesky"
  return c / d, info


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("--terms", required=True)
  ap.add_argument("--gamma", type=float, default=0.0005)
  ap.add_argument("--emax", type=float, default=15000.0)
  ap.add_argument("--ridge", type=float, default=0.0)
  ap.add_argument("--cutoff", type=float, default=0.0)
  ap.add_argument("--test-frac", type=float, default=0.1)
  ap.add_argument("--save", default="")
  a = ap.parse_args()

  crd, Eh, idx = cp.read_energies(f"{EX}/CH3OH-CBS.energies")
  Erel = (Eh - Eh.min()) * cp.HARTREE
  P = np.array(read_terms(a.terms))
  t0 = time.time()
  xi = cp.internal_xi(crd, NL0, "chart")
  F = cp.design_matrix(xi, P)
  print(f"{a.terms}: {F.shape[1]} terms, design matrix {time.time()-t0:.1f}s")
  rng = np.random.default_rng(123)
  test = rng.random(len(Erel)) < a.test_frac
  w = ps_weights(Erel, a.gamma, a.emax)
  wtr = w * (~test)
  t0 = time.time()
  c, info = solve(F, Erel, wtr, ridge=a.ridge, cutoff=a.cutoff)
  res = Erel - F @ c
  print(f"solve {time.time()-t0:.1f}s ({info}); max|c| = {np.abs(c).max():.3e}")
  hdr = " ".join(f"{b:>8.0f}" for b in BANDS)
  print(f"{'':12s}{hdr}")
  print("train rms   " + " ".join(f"{x:8.3f}" for x in band_rms(res, Erel, ~test)))
  print("test  rms   " + " ".join(f"{x:8.3f}" for x in band_rms(res, Erel, test)))
  if a.save:
    np.savez(a.save, P=P, c=c, nl=np.array([NL0[k] for k in NL0]))


if __name__ == "__main__":
  main()
