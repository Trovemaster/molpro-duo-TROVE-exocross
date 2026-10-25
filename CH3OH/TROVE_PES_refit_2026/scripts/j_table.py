#!/usr/bin/env python3
"""j_table.py <rundir> [jmax]: block sizes / symmetry-error check for J=0..jmax runs.

Block sizes are taken from TROVE's "Size of the symmetrized hamiltonian" lines;
"levels" is the number of eigenvalues printed (below the DIAGONALIZER uplimit).
"""
import os
import re
import sys
from trove_energies import parse, count_symmetry_errors

run = sys.argv[1]
jmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
print("| J | A1 block | A2 block | E block | levels printed (A1/A2/E) | symmetry-error lines | End of TROVE |")
print("| :-: | :-: | :-: | :-: | :-: | :-: | :-: |")
for J in range(0, jmax + 1):
  f = os.path.join(run, f"J{J}.out")
  if not os.path.exists(f):
    continue
  size = {}
  done = False
  for l in open(f, errors="replace"):
    m = re.search(r"Size of the symmetrized hamiltonian =\s*(\d+)\s+Symmetry =\s*(\S+)", l)
    if m:
      size[m.group(2)] = int(m.group(1))
    if "End of TROVE" in l:
      done = True
  e = parse(f)
  n = {s: len(e.get(s, [])) for s in ("A1", "A2", "E")}
  print(f"| {J} | {size.get('A1', '–')} | {size.get('A2', '–')} | {size.get('E', '–')} | "
        f"{n['A1']}/{n['A2']}/{n['E']} | {count_symmetry_errors(f)} | {'yes' if done else 'no'} |")
