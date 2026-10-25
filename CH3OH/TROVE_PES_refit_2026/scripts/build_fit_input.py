#!/usr/bin/env python3
"""build_fit_input.py: create CH3OH_fitpes.x / CH3OH_v10 style input files.

Energies are written in cm-1 relative to the lowest ab initio energy
(conv_energy = 1).  Linear terms come from a term list (gen_terms.py) and/or
an existing parameter file (.res/.inp/.trove) used for initial values.
"""
import argparse
import numpy as np
import ch3oh_pes as cp

EX = "/home/user/ExoMol/CH3OH/pes_refit/exomol"
NAMES = ("rCOe", "rOHe", "rCH", "aHOC", "aHCO", "a1", "a2", "b1")


def read_any_params(fname):
  """Parse nonlinear + f-terms from .res/.inp (CH3OH_v10) or .trove files."""
  nl, terms = {}, []
  with open(fname) as f:
    for line in f:
      s = line.split()
      if not s:
        continue
      if s[0].startswith("---") and terms:
        break
      if s[0].lower() == "end" and terms:
        break
      if s[0] in NAMES:
        nl[s[0]] = float(s[-1])
      elif s[0] in ("f", "h") and len(s) >= 15:
        p = tuple(int(x) for x in s[1:13])
        terms.append((p, float(s[14])))
  return nl, terms


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("--terms", help="term list from gen_terms.py")
  ap.add_argument("--init", help="parameter file for initial values")
  ap.add_argument("--data", default=f"{EX}/CH3OH-CBS.energies")
  ap.add_argument("--emax-data", type=float, default=1e9,
                  help="drop points above this E_rel (cm-1)")
  ap.add_argument("--itmax", type=int, default=20)
  ap.add_argument("--factor-out", type=float, default=0.0)
  ap.add_argument("--factor-back", type=float, default=1e-5)
  ap.add_argument("--robust", type=int, default=0)
  ap.add_argument("--tol", type=float, default=-1.0)
  ap.add_argument("--gamma", type=float, default=0.0005)
  ap.add_argument("--ps-emax", type=float, default=15000.0)
  ap.add_argument("--nonlin-ivar", type=int, default=1)
  ap.add_argument("--term-ivar", type=int, default=1)
  ap.add_argument("--ivar-high", type=int, default=0,
                  help="if >0: terms with rigid degree >= this get ivar 1, others ivar 2")
  ap.add_argument("--ridge", type=float, default=0.0)
  ap.add_argument("--tau", default="chart")
  ap.add_argument("--trove-out", default="ch3oh_poten.trove")
  ap.add_argument("--title", default="CH3OH PEF refit to CH3OH-CBS.energies")
  ap.add_argument("--ps-in-python", action="store_true",
                  help="write Partridge-Schwenke weights into the data lines and disable PS in Fortran")
  ap.add_argument("--ps-floor", type=float, default=0.0,
                  help="lower bound of the Partridge-Schwenke weights (with --ps-in-python)")
  ap.add_argument("--pseudo", nargs="*", default=[],
                  help="files with pseudo-points (grid format: 12 coords, E_rel, weight, index)")
  ap.add_argument("--out", required=True)
  a = ap.parse_args()

  nl0 = {"rCOe": 1.42167424, "rOHe": 0.95631449, "rCH": 1.09061196,
         "aHOC": 108.10245297, "aHCO": 110.28972802, "a1": 2.17023887,
         "a2": 2.55392176, "b1": 1.85744195}
  init = {}
  if a.init:
    nli, ti = read_any_params(a.init)
    nl0.update(nli)
    init = {p: v for p, v in ti}
  if a.terms:
    from gen_terms import read_terms
    plist = read_terms(a.terms)
  else:
    plist = list(init.keys())
  zero = tuple([0] * 12)
  if zero not in plist:
    plist = [zero] + plist
  else:
    plist = [zero] + [p for p in plist if p != zero]

  crd, Eh, idx = cp.read_energies(a.data)
  Erel = (Eh - Eh.min()) * cp.HARTREE
  keep = Erel <= a.emax_data
  wdat = np.ones(len(Erel))
  ps_gamma, ps_emax = a.gamma, a.ps_emax
  if a.ps_in_python and a.ps_emax > 0:
    wdat = (np.tanh(-a.gamma * (Erel - a.ps_emax)) + 1.000200020002) / 2.000200020002
    wdat = np.maximum(wdat, a.ps_floor)
    ps_emax = -1.0
  pseudo = []
  for pf in a.pseudo:
    for l in open(pf):
      s = l.split()
      if len(s) >= 14:
        pseudo.append(l.rstrip("\n"))
  with open(a.out, "w") as f:
    f.write("*" * 73 + "\n")
    f.write(f"*  {a.title[:76]}\n*  {len(plist)} linear terms, "
            f"{keep.sum()} points (E_rel <= {a.emax_data:g} cm-1) + {len(pseudo)} pseudo-points\n"
            "*  2026\n*  CH3OH\n")
    f.write("*" * 73 + "\n" + "*" * 73 + "\n*" + " " * 20 + "PARAMETERS OF THE FIT\n" + "*" * 73 + "\n")
    f.write(" 12          <= Number of degrees of freedom\n")
    f.write(" 5           <= Number of bond lengths\n")
    f.write(" 7           <= Number of bond angles\n")
    f.write(f" {a.itmax:<11d} <= Number of iterations in the fit (0-no fit)\n")
    f.write("  1.d-8   1.d-5   <= stab_best, stadev_best\n")
    f.write("  1.0000000000    <= conv_bond\n")
    f.write("  1.0000000000    <= conv_ener (energies in cm-1 relative to Emin)\n")
    f.write(f" {a.factor_out:<12.6f} <= Parameter is out if FACTOR_OUT*|PARAM(I)|<|StErr(I)|\n")
    f.write(f" {a.factor_back:<12.6g} <= Parameter is back if (Delta StDev)/StDev > FACTOR_back\n")
    f.write(f" {a.robust:<12d} <= Watson Robust fit on (1), off (0), or stand. from rob (2)\n")
    f.write(f" {a.tol:<12.4g} <= tol: >0 eigen-cutoff pseudo-inverse, <=0 Cholesky\n")
    f.write(" 1           <= fit_scale\n")
    f.write(f" {ps_gamma:<12.6g} <= Partridge and Schwenke's gamma\n")
    f.write(f" {ps_emax:<12.2f} <= Partridge and Schwenke's Enermax"
            + (" (PS weights written in the data lines)" if a.ps_in_python else "") + "\n")
    if a.ridge > 0:
      f.write(f"ridge {a.ridge:.6e}\n")
    f.write(f"tau {a.tau}\n")
    f.write(f"trove {a.trove_out}\n")
    f.write("*" * 73 + "\n*" + " " * 20 + "POTENTIAL PARAMETERS\n")
    f.write("**********|****|*****|*************|" + "*" * 38 + "\n")
    f.write(" 8                                <= Number of non-linear parameters\n")
    for n in NAMES:
      f.write(f"{n:<10s} {a.nonlin_ivar:4d}   {nl0[n]:20.10f}\n")
    for p in plist:
      v = init.get(p, 0.0)
      if p == zero:
        v = 0.0   # energies are relative to Emin: constant refitted from 0
      deg = sum(p[:11])
      if a.ivar_high > 0:
        iv = 1 if deg >= a.ivar_high else 2
      else:
        iv = a.term_ivar
      if p == zero:
        iv = max(iv, 2) if a.ivar_high > 0 else a.term_ivar
      f.write(" f " + " ".join(f"{x:3d}" for x in p) + f" {iv:3d}  {v:24.14e}\n")
    f.write("end\n")
    for c, e, i, w in zip(crd[keep], Erel[keep], idx[keep], wdat[keep]):
      f.write(" ".join(f"{x:.8f}" for x in c) + f" {e:16.6f} {w:.8e} {i:6d}\n")
    for l in pseudo:
      f.write(l + "\n")
    f.write("-1.0 " * 14 + "\n")
  print(f"wrote {a.out}: {len(plist)} terms, {keep.sum()} points, {len(pseudo)} pseudo-points")


if __name__ == "__main__":
  main()
