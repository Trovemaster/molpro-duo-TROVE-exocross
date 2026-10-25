#!/usr/bin/env python3
"""make_j0contr_inputs.py <rundir> <erange> [jmax]

TROVE J=0-contraction ("model j=0") inputs for a run directory that already
contains a J=0 run (J0.inp/J0.out with EIGENFUNC save, matelem save split):

  J0conv.inp : J=0, model j=0, ERANGE <erange>: CONTRACT save, matelem convert split,
               EIGENFUNC read  -> matrix elements in the J=0 eigenbasis
  J<n>c.inp  : J=n, model j=0, ERANGE <erange>: CONTRACT read, matelem read split,
               EIGENFUNC save

The J=0 eigenstates with E - ZPE <= erange form the vibrational basis for J > 0.
"""
import os
import re
import sys

run, erange = sys.argv[1], float(sys.argv[2])
jmax = int(sys.argv[3]) if len(sys.argv) > 3 else 3
base = open(os.path.join(run, "J0.inp")).read()


def contraction(text):
  # add model j=0 and ERANGE (after ENERCUT, which also sets erange(2)) to CONTRACTION
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


def checkpoint(text, contract, matelem, eigen):
  t = re.sub(r"(?mi)^BASIS_SET\s+\w+", "BASIS_SET   read", text)
  t = re.sub(r"(?mi)^CONTRACT\s+\w+", f"CONTRACT    {contract}", t)
  t = re.sub(r"(?mi)^matelem\s+\w+\s+split", f"matelem     {matelem} split", t)
  t = re.sub(r"(?mi)^EIGENFUNC\s+\w+", f"EIGENFUNC   {eigen}", t)
  return t


t0 = checkpoint(contraction(base), "save", "convert", "read")
t0 = re.sub(r"Jrot\s+\d+", "Jrot  0", t0)
open(os.path.join(run, "J0conv.inp"), "w").write(t0)
for J in range(1, jmax + 1):
  t = checkpoint(contraction(base), "read", "read", "save")
  t = re.sub(r"Jrot\s+\d+", f"Jrot  {J}", t)
  open(os.path.join(run, f"J{J}c.inp"), "w").write(t)
print(f"wrote {run}/J0conv.inp and J1c..J{jmax}c.inp (erange {erange})")
