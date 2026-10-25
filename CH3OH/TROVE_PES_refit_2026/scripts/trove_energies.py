#!/usr/bin/env python3
"""trove_energies.py: parse TROVE 'Variational solution' energies."""
import re
import sys


def parse(fn):
  """Return dict sym -> list of (index, energy) for the last J block."""
  out = {}
  cur = None
  with open(fn) as f:
    for l in f:
      if "Variational solution - irreducible representation" in l:
        cur = "pending"
        continue
      if cur == "pending" and l.strip().startswith("Gamma"):
        cur = "on"
        continue
      if cur == "on":
        s = l.split()
        if len(s) >= 3 and re.match(r"^(A1|A2|E)$", s[0]) and s[1].isdigit():
          out.setdefault(s[0], []).append((int(s[1]), float(s[2])))
        elif not l.strip():
          cur = None
  return out


def count_symmetry_errors(fn):
  n = 0
  with open(fn) as f:
    for l in f:
      ll = l.lower()
      if "different symmetr" in ll or ("symmetr" in ll and "error" in ll):
        n += 1
  return n


if __name__ == "__main__":
  a = parse(sys.argv[1])
  if len(sys.argv) > 2:
    b = parse(sys.argv[2])
    for s in a:
      n = min(len(a[s]), len(b.get(s, [])))
      d = max(abs(x[1] - y[1]) for x, y in zip(a[s][:n], b[s][:n])) if n else float("nan")
      print(f"{s}: {len(a[s])} vs {len(b.get(s, []))} levels; max|dE| (first {n}) = {d:.3e}")
  else:
    for s in a:
      print(s, len(a[s]), [round(e, 4) for _, e in a[s][:10]])
  print("symmetry-error lines:", count_symmetry_errors(sys.argv[1]))
