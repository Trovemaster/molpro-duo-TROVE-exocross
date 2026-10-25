#!/usr/bin/env python3
"""levels_table.py <J0.out> [Emax] [--md]: J=0 levels with largest-contribution labels.

Local-mode quanta (TROVE order): v1 CO, v2 OH, v3-v5 CH, v6 COH, v7-v9 OCH,
v10-v11 S_a/S_b, v12 torsional function index (0 = A1, 1/2 = E of the
ground torsional triplet; higher triplets for larger torsional bases).
"""
import re
import sys

NAMES = ["CO", "OH", "CH", "CH", "CH", "COH", "OCH", "OCH", "OCH", "S", "S"]


def parse(fn):
  out, cur = [], None
  with open(fn) as f:
    for l in f:
      if "Variational solution - irreducible representation" in l:
        cur = "p"
        continue
      if cur == "p" and l.strip().startswith("Gamma"):
        cur = "on"
        continue
      if cur == "on":
        s = l.split()
        if len(s) >= 3 and re.match(r"^(A1|A2|E)$", s[0]) and s[1].isdigit():
          m = re.findall(r"\(([^()]*)\)", l)
          q = tuple(int(x) for x in m[1].split(";")[1].split())
          tail = l.split(")")[2].split("(")[0].split()
          c = float(tail[0]) if tail else float("nan")
          out.append((s[0], int(s[1]), float(s[2]), q, c))
        elif not l.strip():
          cur = None
  return out


def label(q):
  parts = {}
  for i, v in enumerate(q[:11]):
    if v:
      parts[NAMES[i]] = parts.get(NAMES[i], 0) + v
  vib = " + ".join(f"{v if v > 1 else ''}{n}" for n, v in parts.items()) or "ground"
  return vib


if __name__ == "__main__":
  fn = sys.argv[1]
  emax = float(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else 3800.0
  md = "--md" in sys.argv
  lv = sorted(parse(fn), key=lambda t: t[2])
  if md:
    print("| E / cm⁻¹ | Γ | largest contribution (v1 … v12) | |c|² | label |")
    print("| --: | :-: | :-: | --: | --- |")
  for s, i, e, q, c in lv:
    if e > emax:
      break
    qs = " ".join(str(x) for x in q)
    if md:
      print(f"| {e:.2f} | {s} | {qs} | {c:.2f} | {label(q)} |")
    else:
      print(f"{e:10.2f} {s:3s} {qs:32s} {c:5.2f}  {label(q)}")
