#!/usr/bin/env python3
"""make_prod_inputs.py <src J0.inp> <dst dir> --pmax P [--vmax 4] [--tors 8] [--wtors 0.34]
                       [--enercut E] [--uplimit U]

Production-type CH3OH TROVE J=0 input (new PES, exact KEO, R-ALPHA132-THETA-S-TAU, C3v(M)-2):
  * all vibrational classes v <= vmax (CO, OH, CH3 stretches, COH, OCH x3, S_a/S_b);
  * torsion: Fourier, range 0,<tors> (tors+1 functions, i.e. (tors+1)/3 torsional triplets);
  * polyad weights: CO 1, OH 2, CH 2, COH 1, OCH 1, S 1, torsion <wtors>
    (wtors=0.34: one polyad unit per torsional triplet; 0.1: torsion free);
  * PRIMITIVES/CONTRACTION Npolyads = P, CONTRACTION enercut = E (zero-order energy cut);
  * DIAGONALIZER uplimit U.
"""
import argparse
import os
import re
import shutil

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("dst")
ap.add_argument("--pmax", type=int, required=True)
ap.add_argument("--vmax", type=int, default=4)
ap.add_argument("--tors", type=int, default=8)
ap.add_argument("--wtors", type=float, default=0.34)
ap.add_argument("--enercut", type=float, default=80000.0)
ap.add_argument("--uplimit", type=float, default=12000.0)
a = ap.parse_args()
RESC = {1: 1.0, 2: 2.0, 3: 2.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: a.wtors}

t = open(a.src).read()
out, sect = [], None
for line in t.split("\n"):
  s = line.strip()
  u = s.upper()
  if u in ("PRIMITIVES", "CONTRACTION", "BASIS", "DIAGONALIZER"):
    sect = u
  elif sect and u == "END":
    sect = None
  if sect in ("PRIMITIVES", "CONTRACTION") and re.match(r"^NPOLYADS\b", u):
    line = re.sub(r"(?i)(npolyads\s+)\d+", rf"\g<1>{a.pmax}", line)
  if sect == "CONTRACTION" and re.match(r"^ENERCUT\b", u):
    line = re.sub(r"(?i)(enercut\s+)[0-9.]+", rf"\g<1>{a.enercut:.1f}", line)
  if sect == "DIAGONALIZER" and re.match(r"^UPLIMIT\b", u):
    line = re.sub(r"(?i)(uplimit\s+)[0-9.]+", rf"\g<1>{a.uplimit:.1f}", line)
  if sect == "BASIS" and re.match(r"^\d+\s*,", s) and not s.startswith("0"):
    cls = int(s.split(",")[0])
    rmax = a.tors if cls == 7 else a.vmax
    line = re.sub(r"range\s+0\s*,\s*\d+", f"range 0,{rmax}", line)
    line = re.sub(r"resc\s+[0-9.]+", f"resc {RESC[cls]:.2f}", line)
  out.append(line)
os.makedirs(a.dst, exist_ok=True)
open(os.path.join(a.dst, "J0.inp"), "w").write("\n".join(out))
srcdir = os.path.dirname(os.path.abspath(a.src))
for fn in ("kinetic.chk", "combination.chk", "potential.chk", "potential.chk.trove_coords"):
  shutil.copy(os.path.join(srcdir, fn), os.path.join(a.dst, fn))
print(f"wrote {a.dst}/J0.inp: vmax={a.vmax}, tors=0..{a.tors} (w={a.wtors}), Pmax={a.pmax}, "
      f"enercut={a.enercut}, uplimit={a.uplimit}")
