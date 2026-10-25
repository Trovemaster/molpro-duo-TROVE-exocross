#!/usr/bin/env python3
"""summarise_trove.py <rundir>: markdown tables for J=0..3 TROVE runs."""
import os
import sys
from trove_energies import parse, count_symmetry_errors
from compare_trove_j0 import parse_assign

run = sys.argv[1]
old = "/home/user/ExoMol/CH3OH/nb"
print("| J | A1 block | A2 block | E block | total | symmetry-error lines |")
print("| :-: | :-: | :-: | :-: | :-: | :-: |")
for J in range(0, 4):
  f = os.path.join(run, f"J{J}.out")
  if not os.path.exists(f):
    continue
  e = parse(f)
  n = {s: len(e.get(s, [])) for s in ("A1", "A2", "E")}
  print(f"| {J} | {n['A1']} | {n['A2']} | {n['E']} | {sum(n.values())} | {count_symmetry_errors(f)} |")
print()
# J=0: new vs old (assignment matched)
a = parse_assign(os.path.join(run, "J0.out"))
b = parse_assign(os.path.join(old, "J0.out"))
bmap = {}
for s, i, e, q in b:
  bmap.setdefault((s, q), []).append(e)
print("| Γ | (v1 … v12) | E new / cm⁻¹ | E old / cm⁻¹ | new − old |")
print("| :-: | :-: | --: | --: | --: |")
used = {}
for s, i, e, q in sorted(a, key=lambda t: t[2])[:28]:
  k = (s, q); j = used.get(k, 0); used[k] = j + 1
  lst = bmap.get(k, [])
  e2 = lst[j] if j < len(lst) else None
  qs = " ".join(str(x) for x in q)
  if e2 is None:
    print(f"| {s} | {qs} | {e:.3f} | – | – |")
  else:
    print(f"| {s} | {qs} | {e:.3f} | {e2:.3f} | {e - e2:+.3f} |")
print()
# J=1..3 lowest levels (ground vibrational state region < 120 cm-1)
for J in (1, 2, 3):
  f = os.path.join(run, f"J{J}.out")
  if not os.path.exists(f):
    continue
  e = parse(f)
  fo = os.path.join(old, f"J{J}.out")
  eo = parse(fo) if os.path.exists(fo) else {}
  for s in ("A1", "A2", "E"):
    lv = [x for _, x in e.get(s, []) if x < 120.0]
    lo = [x for _, x in eo.get(s, []) if x < 120.0]
    print(f"J={J} {s}: new " + ", ".join(f"{x:.3f}" for x in lv) +
          ("  | old " + ", ".join(f"{x:.3f}" for x in lo) if lo else ""))
