"""A single A-only constrained Duo proposal with independent reassessment.

This diagnostic never promotes a model. It freezes X and all other objects,
uses fixed supplied point weights, and checks a fresh native spectrum rather
than trusting the residual table preceding a proposed parameter update.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import barrier
import coupled


def strict_json(value):
    if isinstance(value, dict):
        return {k: strict_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [strict_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def freeze_other_objects(text):
    edits = []
    for key, f in coupled.fields(text).items():
        if key == ('POTEN', 2, 2):
            continue
        for (_, end), fitted in zip(f.spans, f.fitted):
            if fitted:
                stop = text.index('\n', end)
                edits.append((end, stop))
    for end, stop in sorted(edits, reverse=True):
        text = text[:end]+re.sub(r'(?i)\bfit\b', '', text[end:stop])+text[stop:]
    assert all(not any(f.fitted) for k, f in coupled.fields(text).items() if k != ('POTEN', 2, 2))
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--duo', type=Path, required=True)
    ap.add_argument('--scale', type=float, default=.1)
    ap.add_argument('--jmax', type=int, default=98)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    template = freeze_other_objects(a.input.read_text())
    block = re.search(r'(?ims)^abinitio\s+poten\s+2\b.*?^end\s*$', coupled.uncomment(template))[0]
    if re.search(r'(?im)^\s*weighting\b', block):
        raise ValueError('This diagnostic requires explicit fixed reference point weights')
    values = re.search(r'(?ims)^values\s*\n(.*?)^end\s*$', block)[1]
    ref = np.array([list(map(float, line.split())) for line in values.splitlines() if line.strip()])
    if ref.shape[1] != 3:
        raise ValueError('Expected r, energy and weight columns')
    seed = a.output/'seed.inp'
    seed.write_text(template)
    native = [sys.executable, str(Path(__file__).with_name('run_native.py'))]
    common = ['--duo', str(a.duo), '--jmax', str(a.jmax), '--timeout', '180']
    for filename, name, extra in [(seed, 'baseline', ['--precise']),
                                 (seed, 'proposal', ['--iterations', '1', '--scale', str(a.scale), '--robust', '0'])]:
        result = subprocess.run(native+[str(filename), str(a.output/name)]+common+extra,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        (a.output/f'{name}_console.txt').write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(f'Native run failed: {name}; see console and manifest')
    checkpoints = coupled.checkpoints((a.output/'proposal/duo.out').read_text(errors='replace'))
    if not checkpoints:
        raise ValueError('No complete native proposal')
    candidate = coupled.transfer(template, checkpoints[0])
    (a.output/'candidate.inp').write_text(candidate)
    with np.errstate(over='ignore', invalid='ignore'):
        shape = barrier.audit(candidate)
    result = subprocess.run(native+[str(a.output/'candidate.inp'), str(a.output/'evaluation')]+common+['--precise'],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (a.output/'evaluation_console.txt').write_text(result.stdout)
    if result.returncode:
        raise RuntimeError('Fresh candidate evaluation failed')
    old_fields, new_fields = coupled.fields(template), coupled.fields(candidate)
    for key in old_fields:
        if key != ('POTEN', 2, 2):
            assert old_fields[key].values == new_fields[key].values
    report = dict(status='diagnostic, not selected', fit_scale=a.scale, robust=0,
                  frozen_objects='All except A-state PEC', complete_proposals=len(checkpoints),
                  evaluated_proposal_index_zero_based=0, shape_audit=shape)
    for label, text in [('baseline', template), ('evaluation', candidate)]:
        stats = coupled.save_assessment(a.output/label, a.jmax)
        residual = ref[:, 1]-barrier.potential(ref[:, 0], coupled.fields(text)['POTEN', 2, 2])
        report[label] = dict(experiment=stats,
            reference_weighted_rms_cm=float(np.sqrt(np.sum(ref[:, 2]*residual**2)/sum(ref[:, 2]))))
    report = strict_json(report)
    (a.output/'trial_result.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(shape_ok=shape['ok'],
        baseline_A_rms_cm=report['baseline']['experiment']['by_state']['2']['weighted_rms_cm'],
        trial_A_rms_cm=report['evaluation']['experiment']['by_state']['2']['weighted_rms_cm'],
        reference_before=report['baseline']['reference_weighted_rms_cm'],
        reference_after=report['evaluation']['reference_weighted_rms_cm']), indent=2), flush=True)


if __name__ == '__main__':
    main()
