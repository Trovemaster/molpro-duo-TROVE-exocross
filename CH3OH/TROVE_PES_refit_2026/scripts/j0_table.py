#!/usr/bin/env python3
"""j0_table.py <new J0.out> <old J0.out> [n]: lowest J=0 levels per symmetry."""
import sys
from trove_energies import parse

a = parse(sys.argv[1]); b = parse(sys.argv[2])
n = int(sys.argv[3]) if len(sys.argv) > 3 else 12
print("| # | A1 new | A1 old | A2 new | A2 old | E new | E old |")
print("| :-: | --: | --: | --: | --: | --: | --: |")
for i in range(n):
  row = [str(i + 1)]
  for s in ("A1", "A2", "E"):
    for d in (a, b):
      lv = d.get(s, [])
      row.append(f"{lv[i][1]:.3f}" if i < len(lv) else "–")
  print("| " + " | ".join(row) + " |")
