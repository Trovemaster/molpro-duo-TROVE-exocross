#!/usr/bin/env python3
"""j0_three_way.py <new J0.out> <old-consistent J0.out> <old J0.out> [Emax]

Lowest J=0 levels of the new PES compared with the old PES evaluated
consistently (fitted EQUILIBRIUM/SPECPARAM) and as in the original TROVE run
(inconsistent EQUILIBRIUM/SPECPARAM).  Levels are matched by symmetry and
the quanta of the largest basis contribution, ignoring the torsional label
(v12), whose numbering depends on the 1D torsional basis; the n-th level of
a given (symmetry, v1..v11) in one run is paired with the n-th in the other.
"""
import sys
from compare_trove_j0 import parse_assign


def keyed(levels):
  out, cnt = {}, {}
  for s, _, e, q in sorted(levels, key=lambda t: t[2]):
    k = (s, q[:11])
    n = cnt.get(k, 0)
    cnt[k] = n + 1
    out[k + (n,)] = e
  return out


def main():
  new, cons, old = (parse_assign(f) for f in sys.argv[1:4])
  emax = float(sys.argv[4]) if len(sys.argv) > 4 else 3800.0
  kc, ko = keyed(cons), keyed(old)
  print("| Γ | v1 … v11 (largest contribution) | new PES | old PES, consistent | old PES, original run "
        "| new − old cons. | old orig. − old cons. |")
  print("| :-: | :-: | --: | --: | --: | --: | --: |")
  for k, e in sorted(keyed(new).items(), key=lambda kv: kv[1]):
    if e > emax:
      break
    ec, eo = kc.get(k), ko.get(k)
    f = lambda x: f"{x:.3f}" if x is not None else "–"
    d1 = f"{e - ec:+.2f}" if ec is not None else "–"
    d2 = f"{eo - ec:+.2f}" if (ec is not None and eo is not None) else "–"
    print(f"| {k[0]} | {' '.join(str(x) for x in k[1])} | {e:.3f} | {f(ec)} | {f(eo)} | {d1} | {d2} |")


if __name__ == "__main__":
  main()
