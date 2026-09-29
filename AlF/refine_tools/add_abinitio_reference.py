"""Convert a two-column PEC and add it to a numbered-state Duo input.

Only a vertical shift is applied. The default reference zero is the smallest
published grid value, not a fitted/interpolated minimum. No radial stretching,
energy rescaling, extrapolated observations or invented points are introduced.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
from barrier import potential
from coupled import fields, uncomment

# CODATA 2022; energies in the Wells--Lane supplement already use cm-1.
BOHR_ANGSTROM = 0.529177210544
HARTREE_CM = 219474.63136320


def convert(points, target_te, distance_unit, energy_unit, reference_minimum=None):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError('Expected at least three rows and exactly two columns')
    if not np.isfinite(points).all() or np.any(points[:, 0] <= 0):
        raise ValueError('All values must be finite and distances positive')
    if np.any(np.diff(points[:, 0]) <= 0):
        raise ValueError('Distances must be strictly increasing without duplicates')
    r = points[:, 0] * {'bohr': BOHR_ANGSTROM, 'angstrom': 1.0}[distance_unit]
    e = points[:, 1] * {'hartree': HARTREE_CM, 'cm-1': 1.0}[energy_unit]
    minimum = float(e.min()) if reference_minimum is None else reference_minimum
    shift = target_te - minimum
    return r, e + shift, minimum, shift


def model_minimum(field):
    # Independent of the vibrational ZPE and auxiliary diabatic V0 parameter.
    from scipy.optimize import minimize_scalar
    mesh = np.linspace(1.2, 2.1, 1801)
    i = int(np.argmin(potential(mesh, field)))
    if i in (0, len(mesh)-1):
        raise ValueError('Model inner minimum is outside the supported bracket')
    fit = minimize_scalar(lambda r: float(potential(np.array([r]), field)[0]),
                          bounds=(mesh[i-1], mesh[i+1]), method='bounded',
                          options={'xatol': 1e-13})
    if not fit.success:
        raise ValueError('Could not locate the model inner minimum')
    return float(fit.fun), float(fit.x)


def ps1997_model_weights(text, field, r, alpha, cutoff):
    """Match native Duo's initial-grid PS1997 weights at supplied points.

    Materializing these weights prevents some Duo versions from also assigning
    positive PS1997 weights to automatically generated extrapolation points.
    These weights remain fixed when the analytic PEC parameters change.
    """
    grid = re.search(r'(?ims)^\s*grid\s*$.*?^\s*end\s*$', uncomment(text))[0]
    n = int(re.search(r'(?im)^\s*npoints\s+(\d+)', grid)[1])
    limits = re.search(r'(?im)^\s*range\s+(\S+)\s+(\S+)', grid)
    lo, hi = map(float, limits.groups())
    kind = re.search(r'(?im)^\s*type\s+(\S+)', grid)
    if kind is None or kind[1] != '0':
        raise ValueError('Explicit native PS1997 weights require a uniform type-0 grid')
    mesh = np.linspace(lo, hi, n)
    v = potential(mesh, field)
    # Fortran nint followed by clamping to 1..ngrid, with no +1 in Duo source.
    indices = np.clip(np.floor((r-lo)/((hi-lo)/(n-1)) + .5).astype(int), 1, n) - 1
    return (np.tanh(-alpha*(v[indices] - v.min() - cutoff)) + 1.000020000200002) / 2.000020000200002


def block_text(state, name, lam, mult, r, e, factor, alpha, cutoff, weights=None):
    lines = [f'abinitio poten {state}', f'name "{name}"', f'lambda {lam:g}',
             f'mult {mult:g}', 'type grid', 'units angstrom cm-1',
             f'fit_factor {factor:.12g}']
    if weights is None:
        lines.append(f'Weighting PS1997 {alpha:.12g} {cutoff:.12g}')
    else:
        lines[:0] = [f'(Fixed initial PS1997 weights: alpha={alpha:g}, cutoff={cutoff:g} cm-1)',
                     '(Explicit weights leave native extrapolation points at zero weight.)']
    lines.append('values')
    for i, (x, y) in enumerate(zip(r, e)):
        line = f'{x:.12f} {y:.10f}'
        if weights is not None:
            line += f' {weights[i]:.14e}'
        lines.append(line)
    return '\n'.join(lines + ['end', ''])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model', type=Path)
    ap.add_argument('data', type=Path, help='Two columns: distance, energy; # comments allowed')
    ap.add_argument('output', type=Path, help='New output directory')
    ap.add_argument('--distance-unit', choices=['bohr', 'angstrom'], required=True)
    ap.add_argument('--energy-unit', choices=['hartree', 'cm-1'], required=True)
    ap.add_argument('--state', type=int, default=2)
    ap.add_argument('--te', type=float, help='Target minimum in cm-1; default: actual model minimum')
    ap.add_argument('--reference-minimum', type=float, help='Override reference minimum in cm-1')
    ap.add_argument('--fit-factor', type=float, default=1e2)
    ap.add_argument('--alpha', type=float, default=1e-3)
    ap.add_argument('--cutoff', type=float, default=18000)
    ap.add_argument('--citation', required=True)
    a = ap.parse_args()
    if not all(np.isfinite(x) and x > 0 for x in [a.fit_factor, a.alpha, a.cutoff]):
        ap.error('fit-factor, alpha and cutoff must be positive and finite')
    text = a.model.read_text()
    if re.search(rf'(?im)^\s*abinitio\s+poten\s+{a.state}\b', uncomment(text)):
        ap.error('A reference already exists for this state; refusing to duplicate it')
    field = fields(text)['POTEN', a.state, a.state]
    target, re_model = model_minimum(field)
    if a.te is not None:
        target = a.te
    if not np.isfinite(target) or (a.reference_minimum is not None and not np.isfinite(a.reference_minimum)):
        ap.error('Energy references must be finite')
    r, e, e_min, shift = convert(np.loadtxt(a.data), target, a.distance_unit,
                                a.energy_unit, a.reference_minimum)
    parent = re.search(rf'(?ims)^poten\s+{a.state}\b.*?^end\s*$', uncomment(text))[0]
    lam = float(re.search(r'(?im)^lambda\s+(\S+)', parent)[1])
    mult = float(re.search(r'(?im)^mult\s+(\S+)', parent)[1])
    insert = re.search(r'(?im)^\s*FITTING\s*$', uncomment(text))
    if insert is None:
        ap.error('Expected a FITTING block')
    weights = ps1997_model_weights(text, field, r, a.alpha, a.cutoff)
    a.output.mkdir(parents=True, exist_ok=False)
    comments = (f'(Source: {a.citation})\n'
                f'(A constant {shift:.10f} cm-1 was added to the converted energies.)\n'
                f'(Reference minimum {e_min:.10f}; target Te {target:.10f} cm-1.)\n'
                '(No radial shift or energy rescaling. Original grid points only.)\n')
    for suffix, wt in [('PS1997', None), ('fixed_weights', weights)]:
        block = comments + block_text(a.state, field.name, lam, mult, r, e,
                                      a.fit_factor, a.alpha, a.cutoff, wt)
        (a.output / f'A_reference_{suffix}.inp').write_text(block)
        (a.output / f'AlF_with_reference_{suffix}.inp').write_text(
            text[:insert.start()] + '\n' + block + '\n' + text[insert.start():])
    np.savetxt(a.output/'A_reference_angstrom_cm-1.dat', np.c_[r, e],
               fmt=['%.12f', '%.10f'], header='r_angstrom V_cm-1; '+a.citation)
    result = dict(citation=a.citation, source_data_sha256=hashlib.sha256(a.data.read_bytes()).hexdigest(),
                  source_model_sha256=hashlib.sha256(a.model.read_bytes()).hexdigest(),
                  source_units=[a.distance_unit, a.energy_unit], point_count=len(r),
                  reference_minimum_cm=e_min, reference_minimum_method=(
                      'minimum published grid value' if a.reference_minimum is None else 'explicit override'),
                  target_te_cm=target, model_re_angstrom=re_model, shift_cm=shift,
                  distance_range_angstrom=[float(r[0]), float(r[-1])],
                  last_point_energy_cm=float(e[-1]), fit_factor=a.fit_factor,
                  weighting=dict(type='PS1997', alpha=a.alpha, cutoff_cm=a.cutoff),
                  status='Reference added; model parameters and experimental data unchanged; no refit performed')
    (a.output/'conversion.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
