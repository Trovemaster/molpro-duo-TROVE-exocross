#!/usr/bin/env python3
"""make_bend_inputs.py <src_rundir> <dst_rundir> <Pmax> [jmax]

TROVE CH3OH inputs with the bends (modes 6-9) unfrozen.  Starts from the
J0..J3 inputs of a frozen-bend run (new PES, consistent EQUILIBRIUM/SPECPARAM)
and changes only the basis definition:

  * modes 6 (COH) and 7-9 (OCH): range 0,2 (as for the stretches);
  * polyad weights (resc): CO 1, OH 2, CH 2, COH 1, OCH 1, S_a/S_b 1, torsion 0.1
    (the torsional triplet A1+E of the ground torsional state must not enter
    the polyad, otherwise A and E torsional components would be truncated
    differently);
  * PRIMITIVES/CONTRACTION Npolyads = Pmax.
"""
import os
import re
import shutil
import sys

src, dst, pmax = sys.argv[1], sys.argv[2], int(sys.argv[3])
jmax = int(sys.argv[4]) if len(sys.argv) > 4 else 3
os.makedirs(dst, exist_ok=True)
RESC = {1: 1.0, 2: 2.0, 3: 2.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 0.1}   # by class


def edit(text):
  out, in_basis, in_pc = [], False, None
  for line in text.split("\n"):
    s = line.strip()
    u = s.upper()
    if u == "BASIS":
      in_basis = True
      out.append(line)
      continue
    if in_basis and u == "END":
      in_basis = False
    if u in ("PRIMITIVES", "CONTRACTION"):
      in_pc = u
    elif in_pc and u == "END":
      in_pc = None
    if in_pc and re.match(r"^NPOLYADS\b", u):
      line = re.sub(r"(?i)(npolyads\s+)\d+", rf"\g<1>{pmax}", line)
    if in_basis and re.match(r"^\d+\s*,", s) and not s.startswith("0"):
      cls = int(s.split(",")[0])
      if cls in (4, 5):
        line = re.sub(r"range\s+0\s*,\s*0", "range  0, 2", line)
      line = re.sub(r"resc\s+[0-9.]+", f"resc {RESC[cls]:.1f}", line)
    out.append(line)
  return "\n".join(out)


for J in range(0, jmax + 1):
  with open(os.path.join(src, f"J{J}.inp")) as f:
    t = edit(f.read())
  with open(os.path.join(dst, f"J{J}.inp"), "w") as f:
    f.write(t)
for fn in ("kinetic.chk", "combination.chk", "potential.chk", "potential.chk.trove_coords"):
  shutil.copy(os.path.join(src, fn), os.path.join(dst, fn))
print(f"wrote {dst}/J0..J{jmax}.inp (bends v<=2, Pmax={pmax})")
