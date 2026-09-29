"""Audit the converted Qin grid against native Duo and plot its implications."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from barrier import potential
from coupled import fields, observations, save_assessment
from prepare_qin_reference import extract_pair


def native_points(path):
    return np.array([list(map(float, line.split())) for line in path.read_text().splitlines()
                     if re.match(r'^\s*2\s+[\d.]', line)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('package', type=Path)
    ap.add_argument('--model', required=True, type=Path)
    ap.add_argument('--baseline', required=True, type=Path)
    a = ap.parse_args()
    root = a.package
    inp = root/'inputs'
    data = np.loadtxt(inp/'A_reference_angstrom_cm-1.dat')
    original = (root/'Potential_energy_original.txt').read_text()
    pairs, columns = extract_pair(original)
    raw = np.array([[float(x), float(y)] for x, y, _ in pairs])
    xrows, _ = extract_pair(original, 'X1Sigma+')
    xraw = np.array([[float(x), float(y)] for x, y, _ in xrows])
    conversion = json.loads((inp/'conversion.json').read_text())
    np.testing.assert_allclose(data[:, 0], raw[:, 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(data[:, 1]-raw[:, 1], conversion['shift_cm'], rtol=0, atol=1e-9)
    assert columns == [7, 8] and len(raw) == 66
    # Catch the subtle different grids after X R=1.70 / A R=1.70.
    assert raw[23, 0] == 1.74 and xraw[23, 0] == 1.72
    assert raw[-1, 0] == 9.0 and raw[-1, 1] == -26.42
    baseline_text = a.model.read_text()
    fs = fields(baseline_text)
    te = conversion['target_te_cm']
    d = fs['POTEN', 2, 2].values[15]
    min_i = int(np.argmin(raw[:, 1]))
    outer = np.flatnonzero((raw[:, 0] > 2) & (raw[:, 0] < 4))
    peak_i = outer[np.argmax(raw[outer, 1])]
    spline = CubicSpline(raw[:, 0], raw[:, 1])
    minimum = minimize_scalar(spline, bounds=(1.64, 1.70), method='bounded')
    maximum = minimize_scalar(lambda r: -float(spline(r)), bounds=(2.50, 2.60), method='bounded')
    result = dict(source_columns_one_based=columns, point_count=len(raw),
                  source_sha256=hashlib.sha256((root/'Potential_energy_original.txt').read_bytes()).hexdigest(),
                  sampled_minimum=dict(r_angstrom=float(raw[min_i, 0]), energy_cm=float(raw[min_i, 1])),
                  sampled_barrier=dict(r_angstrom=float(raw[peak_i, 0]), energy_cm=float(raw[peak_i, 1]),
                                       height_above_sampled_minimum_cm=float(raw[peak_i, 1]-raw[min_i, 1]),
                                       height_above_last_point_cm=float(raw[peak_i, 1]-raw[-1, 1])),
                  endpoint=dict(r_angstrom=float(raw[-1, 0]), energy_cm=float(raw[-1, 1]),
                                shifted_energy_cm=float(data[-1, 1]),
                                depth_above_sampled_minimum_cm=float(raw[-1, 1]-raw[min_i, 1]),
                                note='Finite-distance endpoint, not an independently established asymptote'),
                  interpolated_extrema_for_comparison_only=dict(minimum_r=float(minimum.x), minimum_cm=float(minimum.fun),
                                                               maximum_r=float(maximum.x), maximum_cm=float(-maximum.fun)),
                  model_te_cm=te, shift_cm=conversion['shift_cm'],
                  model_fixed_dissociation_cm=d, endpoint_minus_model_dissociation_cm=float(data[-1, 1]-d),
                  source_sampled_A_minus_X_minima_cm=float(raw[min_i, 1]-xraw[:, 1].min()),
                  paper_table2=dict(A_Te_cm=44542.07, A_De_cm=13682.94,
                                    note='Table values disagree with the supplied grid; no reconciliation assumed.'),
                  validation={})
    for variant in ['PS1997', 'fixed_weights']:
        text = (inp/f'AlF_with_reference_{variant}.inp').read_text()
        ff = fields(text)
        assert fs.keys() == ff.keys()
        for key in fs:
            assert fs[key].values == ff[key].values and fs[key].fitted == ff[key].fitted
        assert observations(text) == observations(baseline_text)
        path = root/f'native_{variant}'
        manifest = json.loads((path/'manifest.json').read_text())
        assert manifest['complete_checkpoint'] and manifest['native_diagnostics_present']
        rows = native_points(path/'fit.pot')
        real = np.any(abs(rows[:, 1, None]-data[:, 0]) < 1e-8, axis=1)
        assert real.sum() == len(raw)
        np.testing.assert_allclose(rows[real, 2], data[:, 1], rtol=0, atol=.006)
        assert np.all(rows[real, -1] > 0)
        if variant == 'fixed_weights':
            assert np.all(rows[~real, -1] == 0)
        assert (a.baseline/'rovibronic_energies.dat').read_bytes() == (path/'rovibronic_energies.dat').read_bytes()
        stats = save_assessment(path, 98)
        result['validation'][variant] = dict(native_returncode=manifest['returncode'],
            genuine_points=int(real.sum()), synthetic_points=int((~real).sum()),
            positive_weight_synthetic_points=int(np.count_nonzero(rows[~real, -1])),
            spectrum_identical_to_baseline=True, source_parameters_and_observations_unchanged=True,
            original_weight_experimental_rms_cm=stats['weighted_rms_cm'])

    rr = np.linspace(1.2, 9, 3000)
    model = potential(rr, fs['POTEN', 2, 2])
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5), layout='constrained')
    blue, orange = '#2463a6', '#bd581f'
    ax = axes[0]
    ax.plot(rr, model-te, color=orange, lw=2, label='Current analytic PEC')
    ax.plot(data[:, 0], data[:, 1]-te, '.', color=blue, markersize=5, label='Qin supplied grid, Te aligned')
    ax.plot(data[:, 0], data[:, 1]-te, color=blue, lw=.7, alpha=.7)
    ax.axhline(d-te, color=orange, ls='--', lw=1, label='Current fixed dissociation limit')
    ax.set(xlim=(1.25, 4), ylim=(-500, 23000), xlabel='Internuclear distance / Angstrom',
           ylabel='Energy above current A minimum / cm$^{-1}$', title='A-state well and barrier')
    ax.legend(fontsize=8, loc='upper right')
    ax = axes[1]
    ax.plot(rr, model, color=orange, lw=2, label='Current analytic PEC')
    ax.plot(data[:, 0], data[:, 1], '.-', color=blue, lw=1, markersize=5, label='Qin supplied grid, Te aligned')
    ax.axhline(d, color=orange, ls='--', lw=1)
    ax.annotate(f'At 9 A: {data[-1,1]:,.2f} cm$^{{-1}}$\n{data[-1,1]-d:,.2f} cm$^{{-1}}$ from fixed limit',
                xy=(8.7, data[-1, 1]), xytext=(4.3, 54400), fontsize=9,
                arrowprops=dict(arrowstyle='->', color=blue))
    ax.set(xlim=(3.8, 9.2), ylim=(53200, 56200), xlabel='Internuclear distance / Angstrom',
           ylabel='Energy relative to X minimum / cm$^{-1}$', title='A vertical Te shift leaves a tail mismatch')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Qin 2022 supplied A-state points: conversion and comparison', fontsize=14)
    fig.savefig(root/'reference_comparison.png', dpi=170)
    plt.close(fig)
    (root/'diagnostics.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
