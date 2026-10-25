#!/usr/bin/env python3
"""ch3oh_pes.py: CH3OH PES coordinate/symmetry toolkit.

Reimplements, in vectorised numpy, the coordinate definitions and the
C3v(M) symmetrised expansion of CH3OH_v10_theta12_14-MAPLESym-13-09.f90
(diff_V_tau + ML_symmetry_transformation_XY3_IV), and the TROVE
R-ALPHA132-THETA-S-TAU / potential.chk ("with_modes") conventions.

Ab initio grid columns (Molpro Z-matrix variables, CH3OH-CBS.energies):
  rCO, rOH, rCH1, rCH2, rCH3 [A], aCOH, aOCH1, aOCH2, aOCH3 [deg],
  dHOCH1, dHOCH2, dHOCH3 [deg], E [Eh], index
"""
import numpy as np

PI = np.pi
DEG = PI / 180.0
HARTREE = 2.194746354e5
A_ = 0.5
B_ = 0.5 * np.sqrt(3.0)


# --------------------------------------------------------------------------
# Data and parameter I/O
# --------------------------------------------------------------------------
def read_energies(fname):
  d = np.loadtxt(fname)
  return d[:, :12], d[:, 12], d[:, 13].astype(int)


def read_params(fname):
  """Read a CH3OH_v10 style parameter block (.inp or .res).

  Returns (nonlin: dict name->value, nonlin_list, terms: list of
  (powers[12], ivar, value)). Stops at the first '----' or 'end' line.
  """
  nonlin = {}
  nonlin_order = []
  terms = []
  names = ("rCOe", "rOHe", "rCH", "aHOC", "aHCO", "a1", "a2", "b1")
  with open(fname) as f:
    for line in f:
      s = line.split()
      if not s:
        continue
      if s[0].startswith("---") or s[0].lower() == "end":
        if terms:
          break
        continue
      if s[0] in names and len(s) >= 3:
        nonlin[s[0]] = float(s[2])
        nonlin_order.append(s[0])
      elif s[0] in ("f", "h") and len(s) >= 15:
        p = [int(x) for x in s[1:13]]
        terms.append((p, int(s[13]), float(s[14])))
  return nonlin, nonlin_order, terms


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------
def internal_xi(coords, nl, tau_mode="raw"):
  """Fit coordinates xi(1..12) exactly as in diff_V_tau (CH3OH_v10).

  tau_mode = 'raw'   : tau = (th1+th2+th3-2pi)/3 with raw input dihedrals
                       (as coded in CH3OH_v10).
  tau_mode = 'chart' : tau = th1 + sqrt(2)/3*S_b (mod 2pi), i.e. the value
                       consistent with the TROVE/AB3CD.nb inverse transform
                       t1 = tau - sqrt(2)/3 S_b,
                       t2 = 4pi/3 + tau + S_a/sqrt6 + S_b/(3 sqrt2),
                       t3 = 2pi/3 + tau - S_a/sqrt6 + S_b/(3 sqrt2).
  """
  c = np.atleast_2d(coords)
  rCO, rOH, rH1 = c[:, 0], c[:, 1], c[:, 2]
  rH2, rH3 = c[:, 4], c[:, 3]           # swapped as in CH3OH_v10
  aCOH = c[:, 5] * DEG
  b1_, b2_, b3_ = c[:, 6] * DEG, c[:, 8] * DEG, c[:, 7] * DEG  # swapped
  th1 = np.mod(c[:, 9] * DEG + 2 * PI, 2 * PI)
  th2 = c[:, 11] * DEG                   # swapped
  th3 = c[:, 10] * DEG
  xi = np.empty((c.shape[0], 12))
  xi[:, 0] = 1.0 - np.exp(-nl["a1"] * (rCO - nl["rCOe"]))
  xi[:, 1] = 1.0 - np.exp(-nl["a2"] * (rOH - nl["rOHe"]))
  xi[:, 2] = 1.0 - np.exp(-nl["b1"] * (rH1 - nl["rCH"]))
  xi[:, 3] = 1.0 - np.exp(-nl["b1"] * (rH2 - nl["rCH"]))
  xi[:, 4] = 1.0 - np.exp(-nl["b1"] * (rH3 - nl["rCH"]))
  xi[:, 5] = aCOH - nl["aHOC"] * DEG
  xi[:, 6] = b1_ - nl["aHCO"] * DEG
  xi[:, 7] = b2_ - nl["aHCO"] * DEG
  xi[:, 8] = b3_ - nl["aHCO"] * DEG
  t12 = np.mod(th1 - th2 + 2 * PI, 2 * PI)
  t23 = np.mod(th2 - th3 + 2 * PI, 2 * PI)
  t13 = np.mod(2 * PI - t12 - t23 + 2 * PI, 2 * PI)
  xi[:, 9] = (2.0 * t23 - t12 - t13) / np.sqrt(6.0)
  xi[:, 10] = (t13 - t12) / np.sqrt(2.0)
  if tau_mode == "raw":
    xi[:, 11] = (th1 + th2 + th3 - 2 * PI) / 3.0
  else:
    xi[:, 11] = np.mod(th1 + np.sqrt(2.0) / 3.0 * xi[:, 10], 2 * PI)
  return xi


def chart_to_dihedrals(Sa, Sb, tau):
  """TROVE R-ALPHA132-THETA-S-TAU inverse (returns th1, th2c, th3c)."""
  t1 = tau - np.sqrt(2.0) * Sb / 3.0
  t2 = 4 * PI / 3 + tau + Sa / np.sqrt(6.0) + Sb / (3 * np.sqrt(2.0))
  t3 = 2 * PI / 3 + tau - Sa / np.sqrt(6.0) + Sb / (3 * np.sqrt(2.0))
  return np.mod(t1, 2 * PI), np.mod(t2, 2 * PI), np.mod(t3, 2 * PI)


# --------------------------------------------------------------------------
# C3v(M) operations (ML_symmetry_transformation_XY3_IV of CH3OH_v10)
# --------------------------------------------------------------------------
def sym_ops_xi(xi):
  """Return list of 6 transformed coordinate arrays (npts x 12)."""
  out = []
  e = xi
  # op1 E
  out.append(e.copy())

  def perm(src, pr, pa):
    d = src.copy()
    d[:, 2], d[:, 3], d[:, 4] = src[:, pr[0]], src[:, pr[1]], src[:, pr[2]]
    d[:, 6], d[:, 7], d[:, 8] = src[:, pa[0]], src[:, pa[1]], src[:, pa[2]]
    return d

  # op2: dst3=src5, dst4=src3, dst5=src4 ; S: [[-a,-b],[b,-a]] ; tau+2pi/3
  d = perm(e, (4, 2, 3), (8, 6, 7))
  d[:, 9] = -A_ * e[:, 9] - B_ * e[:, 10]
  d[:, 10] = B_ * e[:, 9] - A_ * e[:, 10]
  d[:, 11] = e[:, 11] + 2 * PI / 3
  out.append(d)
  # op3: dst3=src4, dst4=src5, dst5=src3 ; S: [[-a,b],[-b,-a]] ; tau-2pi/3
  d = perm(e, (3, 4, 2), (7, 8, 6))
  d[:, 9] = -A_ * e[:, 9] + B_ * e[:, 10]
  d[:, 10] = -B_ * e[:, 9] - A_ * e[:, 10]
  d[:, 11] = e[:, 11] - 2 * PI / 3
  out.append(d)
  # op4 (23)*: r2<->r3 ; Sa, -Sb ; -tau
  d = perm(e, (2, 4, 3), (6, 8, 7))
  d[:, 9] = e[:, 9]
  d[:, 10] = -e[:, 10]
  d[:, 11] = -e[:, 11]
  out.append(d)
  # op5: dst3=src5, dst5=src3 ; S: [[-a,-b],[-b,a]] ; -tau-2pi/3
  d = perm(e, (4, 3, 2), (8, 7, 6))
  d[:, 9] = -A_ * e[:, 9] - B_ * e[:, 10]
  d[:, 10] = -B_ * e[:, 9] + A_ * e[:, 10]
  d[:, 11] = -e[:, 11] - 2 * PI / 3
  out.append(d)
  # op6: dst3=src4, dst4=src3 ; S: [[-a,b],[b,a]] ; -tau+2pi/3
  d = perm(e, (3, 2, 4), (7, 6, 8))
  d[:, 9] = -A_ * e[:, 9] + B_ * e[:, 10]
  d[:, 10] = B_ * e[:, 9] + A_ * e[:, 10]
  d[:, 11] = -e[:, 11] + 2 * PI / 3
  out.append(d)
  return out


# --------------------------------------------------------------------------
# Symmetrised basis functions (design matrix) and PES value
# --------------------------------------------------------------------------
def design_matrix(xi, powers, symmetrise=True, norm=6.0):
  """dF[n, i] = (1/6) sum_g prod_m chi_g[m]^p[m] * trig(k*tau_g).

  powers: int array (nterm x 12); powers[:,11] = k (k<0 -> sin(|k| tau)).
  """
  powers = np.asarray(powers, dtype=int)
  npts = xi.shape[0]
  ops = sym_ops_xi(xi) if symmetrise else [xi]
  pmax = int(powers[:, :11].max())
  kmax = int(np.abs(powers[:, 11]).max())
  F = np.zeros((npts, powers.shape[0]))
  for chi in ops:
    # power tables: pw[m][p] -> (npts,)
    pw = np.empty((11, pmax + 1, npts))
    pw[:, 0, :] = 1.0
    for p in range(1, pmax + 1):
      pw[:, p, :] = pw[:, p - 1, :] * chi[:, :11].T
    ck = np.cos(np.outer(np.arange(kmax + 1), chi[:, 11]))
    sk = np.sin(np.outer(np.arange(kmax + 1), chi[:, 11]))
    for i, p in enumerate(powers):
      t = ck[p[11]] if p[11] >= 0 else sk[-p[11]]
      t = t.copy()
      for m in range(11):
        if p[m]:
          t *= pw[m, p[m]]
      F[:, i] += t
  return F / (norm if symmetrise else 1.0)


def eval_pes(coords, nl, terms, tau_mode="raw"):
  xi = internal_xi(coords, nl, tau_mode=tau_mode)
  P = np.array([t[0] for t in terms])
  c = np.array([t[2] for t in terms])
  return design_matrix(xi, P) @ c


# --------------------------------------------------------------------------
# TROVE potential.chk (with_modes) evaluation
# --------------------------------------------------------------------------
def read_potential_chk(fname):
  P, C = [], []
  with open(fname) as f:
    f.readline()
    for line in f:
      s = line.split()
      if int(s[0]) == 987654321:
        break
      C.append(float(s[2]))
      P.append([int(x) for x in s[3:15]])
  return np.array(P, dtype=int), np.array(C)


def trove_chi_pot(coords, eq, beta, tau_mode="chart"):
  """TROVE chi_pot for R-ALPHA132-THETA-S-TAU from grid coordinates.

  eq   : dict re1,re2,re3,alphae1,alphae2 (A, deg)  [EQUILIBRIUM block]
  beta : (beta1, beta2, beta3)                      [SPECPARAM block]
  Uses the same atom association as the fit (H2<->H3 swap of the file).
  """
  nl = {"rCOe": eq["re1"], "rOHe": eq["re2"], "rCH": eq["re3"],
        "aHOC": eq["alphae1"], "aHCO": eq["alphae2"],
        "a1": beta[0], "a2": beta[1], "b1": beta[2]}
  return internal_xi(coords, nl, tau_mode=tau_mode)


def eval_potential_chk(chi, P, C):
  """V = sum_t C_t prod_{m<=11} chi_m^p_m * fourier(p_12, tau)."""
  return design_matrix(chi, P, symmetrise=False) @ C
