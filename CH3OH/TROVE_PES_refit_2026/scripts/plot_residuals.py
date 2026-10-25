#!/usr/bin/env python3
"""plot_residuals.py: obs-calc vs E for a CH3OH_fitpes 00.res file."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else "residuals.png"
rows = []
with open(res) as f:
  for l in f:
    s = l.split()
    if len(s) == 16:
      try:
        rows.append([float(x) for x in s[12:16]])
      except ValueError:
        pass
r = np.array(rows)[:49537]   # ab initio points come first; pseudo-points are appended
E, eps, w = r[:, 0], r[:, 2], r[:, 3]
fig, axs = plt.subplots(1, 2, figsize=(14, 5))
axs[0].plot(E, eps, ".", ms=1.5, alpha=0.5)
axs[0].set_xlabel("E (ab initio) / cm$^{-1}$"); axs[0].set_ylabel("obs - calc / cm$^{-1}$")
axs[0].set_title(f"{len(E)} points (all)")
axs[0].set_ylim(-300, 300)
s = E < 15000
axs[1].plot(E[s], eps[s], ".", ms=1.5, alpha=0.5)
axs[1].set_xlabel("E (ab initio) / cm$^{-1}$"); axs[1].set_ylabel("obs - calc / cm$^{-1}$")
axs[1].set_ylim(-10, 10)
axs[1].set_title("E < 15000 cm$^{-1}$")
for ax in axs:
  ax.axhline(0, color="k", lw=0.5)
fig.tight_layout()
fig.savefig(out, dpi=110)
print("saved", out)
