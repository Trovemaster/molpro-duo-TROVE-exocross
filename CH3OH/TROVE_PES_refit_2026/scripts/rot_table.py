#!/usr/bin/env python3
"""rot_table.py <new rundir> [old rundir] [Ecut]: ground-state rotational levels.

Pairs the levels below Ecut (default 120 cm-1) of each J = 1..3 and symmetry
in energy order and prints a markdown table (new vs old run).
"""
import os
import sys
from trove_energies import parse

new = sys.argv[1]
old = sys.argv[2] if len(sys.argv) > 2 else "/home/user/ExoMol/CH3OH/nb"
ecut = float(sys.argv[3]) if len(sys.argv) > 3 else 120.0
print("| J | Γ | new PES | old, original run | new − old |")
print("| :-: | :-: | --: | --: | --: |")
for J in (1, 2, 3):
  fn, fo = (os.path.join(d, f"J{J}.out") for d in (new, old))
  if not os.path.exists(fn):
    continue
  en = parse(fn)
  eo = parse(fo) if os.path.exists(fo) else {}
  rows = []
  for s in ("A1", "A2", "E"):
    ln = [x for _, x in en.get(s, []) if x < ecut]
    lo = [x for _, x in eo.get(s, []) if x < ecut + 5.0]
    for i, x in enumerate(ln):
      y = lo[i] if i < len(lo) else None
      rows.append((x, s, y))
  for x, s, y in sorted(rows):
    ys = f"{y:.3f}" if y is not None else "–"
    ds = f"{x - y:+.3f}" if y is not None else "–"
    print(f"| {J} | {s} | {x:.3f} | {ys} | {ds} |")
