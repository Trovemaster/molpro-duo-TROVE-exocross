"""Independent uncontracted sinc-DVR check of the uncoupled A(f) channel.

For this input the f parity has no X Sigma+ partner, so this is the exact
single-channel Hamiltonian, up to numerical discretization. The Duo default
diagonal-L2 convention gives J(J+1)-2 for 1Pi when no L2 curve is supplied.
The e parity retains X/A L-uncoupling and is assessed with native Duo instead.
No resonance widths are inferred from the finite-box spectrum.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import time

for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(name, '2')
local_deps = Path(__file__).resolve().parents[1]/'.pecfit_stage'/'plot_deps'
if local_deps.is_dir():
    sys.path.insert(0, str(local_deps))
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from barrier import potential
from coupled import fields

MU = 11.1484731837442
KINETIC = 16.857629194167 / MU
ZPE = 400.147627939163
LIMIT = 55564.02471


def kinetic(r):
    d = np.arange(len(r))[:, None]-np.arange(len(r))[None, :]
    np.fill_diagonal(d, 1)
    matrix = 2.*np.where(d % 2, -1., 1.)/d**2
    np.fill_diagonal(matrix, np.pi**2/3.)
    return matrix * KINETIC/(r[1]-r[0])**2


def effective(r, field, J):
    return potential(np.asarray(r), field) + KINETIC*(J*(J+1)-2)/np.asarray(r)**2


def inner_barrier(field, J):
    r = np.linspace(1.3, 3.8, 5001)
    v = effective(r, field, J)
    indices = np.flatnonzero((v[1:-1] > v[:-2]) & (v[1:-1] > v[2:]))+1
    if not len(indices):
        return None
    i = indices[0]
    fit = minimize_scalar(lambda x: -float(effective(x, field, J)),
                          bounds=(r[i-1], r[i+1]), method='bounded',
                          options={'xatol': 1e-12})
    return float(fit.x), float(-fit.fun)


def spectrum_ceiling(peak, limit=LIMIT):
    """Include levels up to dissociation even when the local barrier lies below it."""
    return max(limit,peak[1])+300 if peak else limit+30000


def spectrum(field, J, r, T, ceiling=None):
    peak = inner_barrier(field, J)
    H = T.copy()
    H.flat[::len(r)+1] += effective(r, field, J)
    stop = ceiling if ceiling is not None else spectrum_ceiling(peak)
    e, psi = eigh(H, subset_by_value=(43000., stop), driver='evr', check_finite=False,
                  overwrite_a=True)
    probability = psi**2
    mean_r = r @ probability
    tail = np.sum(probability[r>r[-1]-1], axis=0)
    inner = np.sum(probability[r<(peak[0] if peak else 2.55)], axis=0)
    result = []
    inner_v = 0
    for idx, energy in enumerate(e):
        localized = inner[idx] > .5
        result.append(dict(J=J, box_index=idx, inner_v=inner_v if localized else '',
                           energy=float(energy-ZPE), below_limit=bool(energy<LIMIT),
                           below_barrier=bool(peak and energy<peak[1]),
                           inner_probability=float(inner[idx]), mean_r=float(mean_r[idx]),
                           tail_probability=float(tail[idx]),
                           flag='b' if tail[idx]<=.01 and mean_r[idx]<=4 else 'u'))
        if localized:
            inner_v += 1
    return result, peak


def write_csv(path, rows):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--jmax', type=int, default=260)
    ap.add_argument('--js', nargs='+', type=int)
    ap.add_argument('--rmax', type=float, default=8.)
    ap.add_argument('--npoints', type=int, default=601)
    a = ap.parse_args()
    fs = fields(a.input.read_text())
    if any(k[0] in ('L2', 'BOBROT', 'SPINORBIT', 'LAMBDAQ') and 2 in k[1:] for k in fs):
        raise ValueError('This isolated A(f) check does not implement additional state-2 operators')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    r = np.linspace(.7, a.rmax, a.npoints)
    T = kinetic(r)
    all_rows, barriers, summary = [], [], []
    begin = time.monotonic()
    for J in (a.js or range(1, a.jmax+1)):
        rows, peak = spectrum(fs['POTEN', 2, 2], J, r, T)
        all_rows.extend(rows)
        barriers.append(dict(J=J, barrier_r=peak[0] if peak else '',
                             barrier_energy=peak[1]-ZPE if peak else ''))
        bound = [x for x in rows if x['below_limit'] and x['inner_probability']>.5]
        trapped = [x for x in rows if x['below_barrier'] and x['inner_probability']>.5]
        summary.append(dict(J=J, inner_bound_count=len(bound),
                            below_barrier_localized_count=len(trapped),
                            outer_below_limit_count=sum(x['below_limit'] and x['inner_probability']<=.5 for x in rows),
                            duo_style_localized_count=sum(x['flag']=='b' for x in rows)))
        if J % 20 == 0:
            print(json.dumps(summary[-1]), flush=True)
    write_csv(out/'spectrum.csv', all_rows)
    write_csv(out/'barriers.csv', barriers)
    write_csv(out/'counts.csv', summary)
    (out/'manifest.json').write_text(json.dumps(dict(input=str(a.input.resolve()),
        rmax=a.rmax, npoints=a.npoints, Js=a.js or list(range(1, a.jmax+1)),
        kinetic_coefficient=KINETIC, zpe=ZPE, dissociation=LIMIT-ZPE,
        elapsed_seconds=time.monotonic()-begin, status='complete'), indent=2))
    print(json.dumps(dict(last_bound_J=max((x['J'] for x in summary if x['inner_bound_count']),default=None),
                          last_localized_J=max((x['J'] for x in summary if x['below_barrier_localized_count']),default=None))))


if __name__ == '__main__':
    main()
