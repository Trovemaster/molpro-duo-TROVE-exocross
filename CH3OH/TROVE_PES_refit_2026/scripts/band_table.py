#!/usr/bin/env python3
"""band_table.py label=J0.out [label=J0.out ...]: J=0 band-origin table (markdown).

States are picked by symmetry, the largest-contribution label (v1..v11, see
levels_table.py) and an energy window; torsional states use the ground
vibrational label.  Experimental values: Shimanouchi/NIST fundamentals and
CH3OH torsional data (approximate).
"""
import sys
from levels_table import parse, label

# (name, exp, symmetry, label test, window, pick)
def vib(q):
  return tuple(q[:11])

TARGETS = [
  ("v_t=0 tunnelling", "9.12", "E", lambda q: sum(q[:11]) == 0, (0, 60), 0),
  ("v_t=1", "≈ 209", "E", lambda q: sum(q[:11]) == 0, (100, 400), 0),
  ("v_t=1", "≈ 294", "A2", lambda q: sum(q[:11]) == 0, (100, 450), 0),
  ("v_t=2", "≈ 353", "A1", lambda q: sum(q[:11]) == 0, (100, 500), 0),
  ("ν8 CO stretch", "1033", "A1", lambda q: label(q) == "CO", (950, 1200), 0),
  ("ν7 CH3 rock a′", "1060", "A1", lambda q: label(q) in ("OCH", "COH"), (1050, 1200), 0),
  ("ν11 CH3 rock a″", "1165", "A2", lambda q: label(q) == "OCH", (1100, 1300), 0),
  ("ν6 OH bend", "1345", "A1", lambda q: label(q) == "COH", (1300, 1420), 0),
  ("ν5 CH3 umbrella", "1455", "A1", lambda q: label(q) == "OCH" and q[6] == 1, (1400, 1650), 0),
  ("ν10 CH3 deform. a″", "1477", "A2", lambda q: label(q) == "S", (1440, 1560), 0),
  ("ν4 CH3 deform. a′", "1477", "A1", lambda q: label(q) == "S", (1440, 1560), 0),
  ("ν3 CH3 s-stretch", "2844", "A1", lambda q: label(q) == "CH", (2800, 2950), 0),
  ("ν9 CH3 a-stretch a″", "2960", "A2", lambda q: label(q) == "CH", (2900, 3000), 0),
  ("ν2 CH3 a-stretch a′", "3000", "A1", lambda q: label(q) == "CH", (2980, 3070), 0),
  ("ν1 OH stretch", "3681", "A1", lambda q: label(q) == "OH", (3600, 3760), 0),
]

runs = []
for arg in sys.argv[1:]:
  name, fn = arg.split("=", 1)
  runs.append((name, sorted(parse(fn), key=lambda t: t[2])))

print("| band | exp. | Γ | " + " | ".join(n for n, _ in runs) + " |")
print("| --- | --: | :-: | " + " | ".join("--:" for _ in runs) + " |")
for name, exp, sym, test, (lo, hi), pick in TARGETS:
  row = [name, exp, sym]
  for _, lv in runs:
    cand = [e for s, i, e, q, c in lv if s == sym and lo <= e <= hi and test(q)]
    row.append(f"{cand[pick]:.2f}" if len(cand) > pick else "–")
  print("| " + " | ".join(row) + " |")
