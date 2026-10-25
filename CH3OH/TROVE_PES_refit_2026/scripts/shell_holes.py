#!/usr/bin/env python3
"""shell_holes.py: hole search inside the "physical shell" V_sep < Ecap.

V_sep(q) = V(q0) + sum_m [V(q0 with q_m) - V(q0)] is the separable estimate
built from the (well determined) 1D cuts of the fit.  Holes are points with
V(q) < Vlow or V(q) < frac*V_sep(q) (for V_sep > Emin_frac).

Search: (i) random sampling in the shell, (ii) projected descent on V that
rejects steps leaving the shell.  Writes the hole geometries (fit coords)
to an .npz file for patching.
"""
import argparse
import numpy as np
import hole_check as hc

DEG = hc.DEG


def vsep(pef, q, q0, V0):
  n = q.shape[0]
  allq = np.empty((11 * n, 12))
  for m in range(11):
    qq = np.tile(q0, (n, 1))
    qq[:, m] = q[:, m]
    qq[:, 11] = q[:, 11]
    allq[m * n:(m + 1) * n] = qq
  v = pef.V_xi(pef.phys_to_xi(allq)).reshape(11, n)
  vt = pef.V_xi(pef.phys_to_xi(np.column_stack([np.tile(q0[:11], (n, 1)), q[:, 11]])))
  # torsional part enters through vt; 1D increments relative to the same tau
  return vt + np.sum(v - vt[None, :], axis=0)


def sample(pef, n, rng, q0, scale=0.45):
  q = np.tile(q0, (n, 1))
  for m in range(11):
    lo, hi = hc.BOX[m]
    q[:, m] = np.clip(q0[m] + rng.normal(0, scale * (hi - lo), n) * rng.random(n), lo, hi)
  q[:, 11] = rng.uniform(0, 2 * np.pi, n)
  return q


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("pef")
  ap.add_argument("--ecap", type=float, default=30000.0)
  ap.add_argument("--nrand", type=int, default=400000)
  ap.add_argument("--nmin", type=int, default=3000)
  ap.add_argument("--steps", type=int, default=150)
  ap.add_argument("--frac", type=float, default=0.3)
  ap.add_argument("--vlow", type=float, default=-10.0)
  ap.add_argument("--out", default="holes.npz")
  ap.add_argument("--seed", type=int, default=7)
  a = ap.parse_args()
  pef = hc.load(a.pef)
  rng = np.random.default_rng(a.seed)
  q0 = pef.ref_q()
  V0 = pef.V_xi(pef.phys_to_xi(q0[None, :]))[0]
  q = sample(pef, a.nrand, rng, q0)
  V = pef.V_xi(pef.phys_to_xi(q))
  Vs = vsep(pef, q, q0, V0)
  shell = Vs < a.ecap
  bad = shell & (Vs > 1500) & ((V < a.vlow) | ((V < a.frac * Vs) & (Vs > 5000)))
  print(f"random: {shell.sum()} of {a.nrand} in shell V_sep<{a.ecap:.0f};"
        f" holes: {bad.sum()} (V<{a.vlow}: {np.sum(shell & (V < a.vlow))});"
        f" min V in shell = {V[shell].min():.1f}", flush=True)
  holes = [q[bad]]
  # projected descent inside the shell
  qs = q[shell][: a.nmin].copy()
  Vq = V[shell][: a.nmin].copy()
  lo = np.array([b[0] for b in hc.BOX]); hi = np.array([b[1] for b in hc.BOX])
  step = np.array([0.02] * 5 + [0.02] * 6 + [0.05])
  h = 1e-3
  for it in range(a.steps):
    G = np.zeros_like(qs)
    allq = np.concatenate([qs + h * np.eye(12)[m] for m in range(12)] +
                          [qs - h * np.eye(12)[m] for m in range(12)])
    vv = pef.V_xi(pef.phys_to_xi(allq)).reshape(24, -1)
    G = ((vv[:12] - vv[12:]) / (2 * h)).T
    gn = np.sqrt(np.sum((G * step) ** 2, axis=1)) + 1e-30
    qn = np.clip(qs - (G * step) / gn[:, None] * step, lo, hi)
    Vn = pef.V_xi(pef.phys_to_xi(qn))
    Vsn = vsep(pef, qn, q0, V0)
    ok = (Vn < Vq) & (Vsn < a.ecap)
    qs[ok] = qn[ok]; Vq[ok] = Vn[ok]
    step *= 1.02 if ok.mean() > 0.5 else 0.85
    if it % 30 == 0 or it == a.steps - 1:
      print(f"  descent it {it:4d}: min V = {Vq.min():12.1f}; below {a.vlow}: {np.sum(Vq < a.vlow)};"
            f" accepted {ok.mean():.2f}", flush=True)
  Vs_end = vsep(pef, qs, q0, V0)
  badm = (Vs_end > 1500) & ((Vq < a.vlow) | ((Vq < a.frac * Vs_end) & (Vs_end > 5000)))
  holes.append(qs[badm])
  print(f"descent: {badm.sum()} of {len(qs)} end points are holes; lowest V = {Vq.min():.1f}")
  if badm.any():
    i = np.argmin(Vq)
    qq = qs[i].copy(); qq[5:] /= DEG
    print("  lowest point (rCO,rOH,r1,r2,r3 A; aCOH,b1,b2,b3,Sa,Sb,tau deg):")
    print("  ", np.array2string(qq, precision=3), f" V_sep = {Vs_end[i]:.0f}")
  H = np.concatenate(holes)
  np.savez(a.out, q=H, V=pef.V_xi(pef.phys_to_xi(H)) if len(H) else np.zeros(0),
           Vsep=vsep(pef, H, q0, V0) if len(H) else np.zeros(0))
  print(f"saved {len(H)} hole geometries to {a.out}")


if __name__ == "__main__":
  main()
