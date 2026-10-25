#!/usr/bin/env python3
"""make_slice_inputs.py <rundir> <erange> <zpe or 0> [jmax]

TROVE inputs for the sliced J=0-contraction workflow, starting from <rundir>/J0.inp
(a J=0 input with BASIS_SET/CONTRACT/matelem save):

  J0s.inp         J=0: matelem save split 0 0 (vibrational part only), EIGENFUNC save
  slice_S.inp     J=0: matelem save split S S (one rotational/Coriolis slice, S=1..12),
                  EIGENFUNC none, BASIS_SET/CONTRACT read  -> matelemS.chk (N x N)
  conv_S.inp      J=0: model j=0, ERANGE <erange>, matelem convert split S S, EIGENFUNC read,
                  CONTRACT save                      -> j0_matelemS.chk (Nroots x Nroots)
  J<n>c.inp       J=n: model j=0, ERANGE <erange>, CONTRACT read, matelem read split, EIGENFUNC save
                  (DIAGONALIZER ZPE <zpe> if zpe > 0)

Slices: 1-9 = G_rot(a,b), 10-12 = G_cor(a) (TROVE 'split' layout for 12 modes).
"""
import os
import re
import sys

run, erange, zpe = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
jmax = int(sys.argv[4]) if len(sys.argv) > 4 else 3
base = open(os.path.join(run, "J0.inp")).read()


def contraction(text):
  out, inside = [], False
  for line in text.split("\n"):
    u = line.strip().upper()
    if u == "CONTRACTION":
      inside = True
    elif inside and u == "END":
      out.append("  model     j=0")
      out.append(f"  erange    {erange:.1f}")
      inside = False
    out.append(line)
  return "\n".join(out)


def checkpoint(text, basis, contract, matelem, eigen):
  t = re.sub(r"(?mi)^BASIS_SET\s+\w+", f"BASIS_SET   {basis}", text)
  t = re.sub(r"(?mi)^CONTRACT\s+\w+", f"CONTRACT    {contract}", t)
  t = re.sub(r"(?mi)^matelem\s+.*$", f"matelem     {matelem}", t)
  t = re.sub(r"(?mi)^EIGENFUNC\s+\w+", f"EIGENFUNC   {eigen}", t)
  return t


def jrot(text, j):
  return re.sub(r"Jrot\s+\d+", f"Jrot  {j}", text)


def zpe_line(text):
  if zpe <= 0:
    return text
  return re.sub(r"(?mi)^(\s*uplimit\s+[0-9.]+\s*)$", rf"\1\n  ZPE     {zpe:.6f}", text, count=1)


w = lambda name, t: open(os.path.join(run, name), "w").write(t)
w("J0s.inp", jrot(checkpoint(base, "save", "save", "save split 0 0", "save"), 0))
for s in range(1, 13):
  w(f"slice_{s}.inp", jrot(checkpoint(base, "read", "read", f"save split {s} {s}", "none"), 0))
  w(f"conv_{s}.inp", jrot(checkpoint(contraction(base), "read", "save", f"convert split {s} {s}", "read"), 0))
for J in range(1, jmax + 1):
  w(f"J{J}c.inp", zpe_line(jrot(checkpoint(contraction(base), "read", "read", "read split", "save"), J)))
print(f"wrote J0s.inp, slice_1..12.inp, conv_1..12.inp, J1c..J{jmax}c.inp in {run}")
