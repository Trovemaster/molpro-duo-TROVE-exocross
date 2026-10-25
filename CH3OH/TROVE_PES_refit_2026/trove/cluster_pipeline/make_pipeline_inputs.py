#!/usr/bin/env python3
"""make_pipeline_inputs.py <base J0.inp> <outdir> <prefix> --erange E [--jmax J] [--zpe Z]

Writes the TROVE inputs for a CH3OH J=0-contraction ("model j=0") pipeline in which all
steps run in one working directory (e.g. consecutive batch jobs on a cluster):

  <prefix>_J0s.inp     J=0: BASIS_SET save, CONTRACT save, matelem save split 0 0,
                       EIGENFUNC save  -> J=0 energies, contr_matelem.chk, eigen_*0_*.chk
  <prefix>_slices.inp  J=0: BASIS_SET read, CONTRACT read, matelem save split 1 12,
                       EIGENFUNC none  -> matelem1..12.chk (N x N each)
  <prefix>_J0conv.inp  J=0, model j=0, ERANGE E: CONTRACT save, matelem convert split,
                       EIGENFUNC read  -> j0contr_*.chk, j0_matelem*.chk
  <prefix>_J<n>c.inp   J=n, model j=0, ERANGE E: CONTRACT read, matelem read split,
                       EIGENFUNC save, DIAGONALIZER ZPE Z  (only if --zpe is given)

The base input must be a J=0 input whose CHECK_POINT block has BASIS_SET, CONTRACT,
matelem ... split and EIGENFUNC lines (e.g. make_prod_inputs.py output).
Slices: 0 = vibrational Hamiltonian, 1-9 = G_rot(a,b), 10-12 = G_cor(a).
ERANGE is placed after ENERCUT (ENERCUT also sets erange(2) in TROVE).
"""
import argparse
import os
import re

ap = argparse.ArgumentParser()
ap.add_argument("base")
ap.add_argument("outdir")
ap.add_argument("prefix")
ap.add_argument("--erange", type=float, required=True)
ap.add_argument("--jmax", type=int, default=3)
ap.add_argument("--zpe", type=float, default=None)
a = ap.parse_args()

base = open(a.base).read()
os.makedirs(a.outdir, exist_ok=True)


def set_checkpoint(text, basis, contract, matelem, eigen):
  t = re.sub(r"(?mi)^BASIS_SET\s+\w+.*$", f"BASIS_SET   {basis}", text)
  t = re.sub(r"(?mi)^CONTRACT\s+\w+.*$", f"CONTRACT    {contract}", t)
  t = re.sub(r"(?mi)^matelem\s+\w+.*$", f"matelem     {matelem}", t)
  t = re.sub(r"(?mi)^EIGENFUNC\s+\w+.*$", f"EIGENFUNC   {eigen}", t)
  for key in ("BASIS_SET", "CONTRACT", "matelem", "EIGENFUNC"):
    if not re.search(rf"(?mi)^{key}\s", t):
      raise SystemExit(f"CHECK_POINT line {key} not found in {a.base}")
  return t


def set_jrot(text, j):
  t, n = re.subn(r"Jrot\s+\d+", f"Jrot  {j}", text)
  if n != 1:
    raise SystemExit(f"expected exactly one 'Jrot' in {a.base}, found {n}")
  return t


def add_j0_model(text):
  out, inside, done = [], False, False
  for line in text.split("\n"):
    u = line.strip().upper()
    if u == "CONTRACTION":
      inside = True
    elif inside and u == "END":
      out.append("  model     j=0")
      out.append(f"  erange    {a.erange:.1f}")
      inside, done = False, True
    out.append(line)
  if not done:
    raise SystemExit("CONTRACTION block not found")
  return "\n".join(out)


def add_zpe(text, zpe):
  out, inside, done = [], False, False
  for line in text.split("\n"):
    u = line.strip().upper()
    if u == "DIAGONALIZER":
      inside = True
    elif inside and u == "END":
      out.append(f"  ZPE     {zpe:.6f}")
      inside, done = False, True
    out.append(line)
  if not done:
    raise SystemExit("DIAGONALIZER block not found")
  return "\n".join(out)


def write(name, text):
  path = os.path.join(a.outdir, f"{a.prefix}_{name}.inp")
  open(path, "w").write(text)
  print("wrote", path)


write("J0s", set_jrot(set_checkpoint(base, "save", "save", "save split 0 0", "save"), 0))
write("slices", set_jrot(set_checkpoint(base, "read", "read", "save split 1 12", "none"), 0))
write("J0conv", set_jrot(set_checkpoint(add_j0_model(base), "read", "save", "convert split",
                                        "read"), 0))
if a.zpe is not None:
  for j in range(1, a.jmax + 1):
    t = set_checkpoint(add_j0_model(base), "read", "read", "read split", "save")
    write(f"J{j}c", add_zpe(set_jrot(t, j), a.zpe))
else:
  print("no --zpe: J>0 inputs not written (run again with the J=0 ZPE)")
