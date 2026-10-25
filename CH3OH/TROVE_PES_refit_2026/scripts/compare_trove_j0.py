#!/usr/bin/env python3
"""compare_trove_j0.py: compare TROVE J=0 band origins of two runs by assignment."""
import re
import sys


def parse_assign(fn):
  out = []
  cur = None
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
          # quanta of the largest contribution: 12 integers after ';' in 2nd ()
          m = re.findall(r"\(([^()]*)\)", l)
          q = None
          if len(m) >= 2 and ";" in m[1]:
            q = tuple(int(x) for x in m[1].split(";")[1].split())
          out.append((s[0], int(s[1]), float(s[2]), q))
        elif not l.strip():
          cur = None
  return out


if __name__ == "__main__":
  a = parse_assign(sys.argv[1])
  b = parse_assign(sys.argv[2])
  nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 25
  bmap = {}
  for s, i, e, q in b:
    bmap.setdefault((s, q), []).append(e)
  allst = sorted(a, key=lambda t: t[2])[:nmax]
  print(f"{'sym':3s} {'quanta (v1..v12)':40s} {'E(run1)':>12s} {'E(run2)':>12s} {'diff':>9s}")
  used = {}
  for s, i, e, q in allst:
    k = (s, q)
    j = used.get(k, 0)
    e2 = bmap.get(k, [None] * (j + 1))
    e2 = e2[j] if j < len(e2) else None
    used[k] = j + 1
    qs = " ".join(str(x) for x in q) if q else "?"
    if e2 is None:
      print(f"{s:3s} {qs:40s} {e:12.4f} {'-':>12s}")
    else:
      print(f"{s:3s} {qs:40s} {e:12.4f} {e2:12.4f} {e - e2:9.3f}")
