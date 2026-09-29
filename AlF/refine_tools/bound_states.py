"""Run Duo's localization diagnostic without fitting; retain native evidence.

The b/u label is a finite-box localization test, not a resonance width or
a test that the energy lies below the dissociation threshold.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from coupled import fields, observations, uncomment


def prepare(text, jmax=40, rmax=8., npoints=1001, vmax=100,
            density=1e-2, mean_r=4., energy_max=70000.):
    # In Duo the energies list and FITTING share a single terminating END.
    pattern = r'(?ims)^\s*FITTING\s*$.*?^\s*energies\s*$.*?^\s*end\s*$'
    matches = list(re.finditer(pattern, uncomment(text)))
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one FITTING/energies block, got {len(matches)}')
    match = matches[0]
    text = text[:match.start()] + '\n' + text[match.end():]
    if re.search(r'(?im)^\s*INTENSITY\s*$', text):
        raise ValueError('Input already contains INTENSITY')
    text = re.sub(r'(?im)^\s*jrot\s+[^\n]*', f'jrot 0 - {jmax}', text, count=1)
    text = re.sub(r'(?im)^\s*npoints\s+[^\n]*', f' npoints {npoints}', text, count=1)
    text = re.sub(r'(?im)^\s*range\s+[^\n]*', f' range 0.7 {rmax:g}', text, count=1)
    text = re.sub(r'(?im)^\s*vmax\s+[^\n]*', f' vmax {vmax} {vmax}', text, count=1)
    # Fit annotations are not needed in an intensity-only calculation.
    text = re.sub(r'(?im)\s+fit[ \t]*$', '', text)
    return text + f'''
INTENSITY
absorption
unbound
states_only
THRESH_DELTA_R 1
THRESH_BOUND {density:g}
THRESH_INTES 1e-80
THRESH_LINE 1e-60
TEMPERATURE 2000
thresh_dipole 1e-9
THRESH_bound_rmax {mean_r:g}
print_bound_density
linelist AlF_refined_J{jmax}_unbound
J, 0, {jmax}
freq-window 0, {energy_max:g}
energy low -0.001, 56000.00, upper -0.00, {energy_max:g}
END
'''


def read_states(path):
    """Parse this build's 14-column non-Lande, non-matelem .states format."""
    rows = []
    with Path(path).open() as f:
        for line_number, line in enumerate(f, 1):
            p = line.split()
            if len(p) != 14 or p[11] not in ('b', 'u'):
                raise ValueError(f'Unexpected .states format at {path}:{line_number}: {line}')
            rows.append(dict(id=int(p[0]), energy=float(p[1]), g=int(p[2]), J=int(p[3]),
                             parity=p[4], ef=p[5], state=p[6], v=int(p[7]),
                             Lambda=int(p[8]), Sigma=float(p[9]), Omega=float(p[10]),
                             flag=p[11], mean_r=float(p[12]), tail_probability=float(p[13])))
    if not rows or len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Empty file or duplicate state identifiers')
    return rows


def summarize(source, rows, zpe, output, upper_name='A1Pi', dissociation=None):
    """AlF summary. Keep threshold tests separate from the energy criterion."""
    f = fields(source)['POTEN', 2, 2]
    if dissociation is not None:
        if not math.isfinite(dissociation):
            raise ValueError('Dissociation energy must be finite')
        limit = dissociation-zpe
    elif f.kind != 'COUPLED-PEC':
        raise ValueError('The AlF summary expects the supplied coupled-pec state 2')
    else:
        n1 = int(f.values[7])+9
        limit = f.values[n1]-zpe
    upper = [dict(r) for r in rows if r['state']==upper_name]
    for r in upper:
        r['below_dissociation'] = r['energy'] < limit
        r['category'] = ('below_limit_localized' if r['flag']=='b' else 'below_limit_diffuse') \
            if r['below_dissociation'] else ('above_limit_localized_candidate' if r['flag']=='b' else 'above_limit_box_state')
    categories = {c: sum(r['category']==c for r in upper) for c in sorted({r['category'] for r in upper})}
    localized = [r for r in upper if r['category']=='below_limit_localized']
    obs = [x for x in observations(source) if x['state']==2 and x['J']>=0 and x['observed_cm']>=0]
    by_v = []
    for v in sorted({r['v'] for r in localized} | {r['v'] for r in obs}):
        subset = [r for r in localized if r['v']==v]
        measured = [r for r in obs if r['v']==v]
        by_v.append(dict(v=v, observed_Jmax=max((r['J'] for r in measured),default=None),
                         localized_below_limit_Jmax=max((r['J'] for r in subset),default=None),
                         count_e=sum(r['ef']=='e' for r in subset), count_f=sum(r['ef']=='f' for r in subset)))
    lookup = {}
    for r in upper:
        lookup.setdefault((r['J'],r['parity'],r['v']),[]).append(r)
    checked=[]
    calculated_jmax=max(r['J'] for r in rows)
    for o in obs:
        matches=lookup.get((o['J'],o['parity'],o['v']),[])
        if o['J']>calculated_jmax:
            continue
        if len(matches)!=1:
            raise ValueError(f'Ambiguous/missing observed assignment: {o}, candidates={len(matches)}')
        c=matches[0]
        checked.append(dict(J=o['J'],parity=o['parity'],v=o['v'],observed=o['observed_cm'],
                            calculated=c['energy'],flag=c['flag'],mean_r=c['mean_r'],
                            tail_probability=c['tail_probability'],below_dissociation=c['below_dissociation']))
    result=dict(dissociation_above_X_ground_cm=limit,zpe_cm=zpe,upper_state=upper_name,
                categories=categories,by_v=by_v,observed_checked=len(checked),
                all_checked_observations_localized_below_limit=all(r['flag']=='b' and r['below_dissociation'] for r in checked),
                Jmax_calculated=max(r['J'] for r in rows),
                J_cutoff_reached_by_localized_below_limit_levels=any(r['J']==calculated_jmax for r in localized),
                counting_convention='Each e/f rovibronic level counted once; excludes M and nuclear-spin degeneracy.',
                warning='b/u is localization only, not an inner-well assignment: outer-well and mixed levels may also be b. High printed v values include box-state ordering. Above-limit candidates are not widths or established resonances.')
    for filename, data in [('A_states.csv',upper),('observed_bound_check.csv',checked),('by_v.csv',by_v)]:
        if data:
            with (Path(output)/filename).open('w',newline='') as handle:
                w=csv.DictWriter(handle,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    (Path(output)/'bound_summary.json').write_text(json.dumps(result,indent=2))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--duo', default=os.environ.get('DUO_EXE', 'duo'))
    ap.add_argument('--jmax', type=int, default=40)
    ap.add_argument('--rmax', type=float, default=8.)
    ap.add_argument('--npoints', type=int, default=1001)
    ap.add_argument('--vmax', type=int, default=100)
    ap.add_argument('--density', type=float, default=1e-2)
    ap.add_argument('--mean-r', type=float, default=4.)
    ap.add_argument('--energy-max', type=float, default=70000.)
    ap.add_argument('--dissociation', type=float,
                    help='Absolute asymptote in cm-1 on the input PEC energy zero; required for GRID A')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=7200.)
    a = ap.parse_args()
    duo = Path(shutil.which(a.duo) or a.duo).resolve()
    if not duo.is_file():
        ap.error('Duo executable not found')
    source = a.input.resolve()
    prepared = prepare(source.read_text(), a.jmax, a.rmax, a.npoints, a.vmax,
                       a.density, a.mean_r, a.energy_max)
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    inp = out / 'input.inp'
    inp.write_text(prepared)
    manifest = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    manifest.update(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(), duo=str(duo),
                    duo_sha256=hashlib.sha256(duo.read_bytes()).hexdigest(), status='running')
    env = os.environ.copy()
    paths = ['C:/Program Files (x86)/Intel/oneAPI/compiler/2025.3/bin',
             'C:/Program Files (x86)/Intel/oneAPI/mkl/2025.3/bin']
    env['PATH'] = os.pathsep.join([p for p in paths if Path(p).is_dir()] + [env['PATH']])
    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        env[name] = str(a.threads)
    begin = time.monotonic()
    mp = out / 'manifest.json'
    mp.write_text(json.dumps(manifest, indent=2))
    with inp.open('rb') as stdin, (out / 'duo.out').open('wb') as stdout:
        proc = subprocess.Popen([str(duo)], stdin=stdin, stdout=stdout,
                                stderr=subprocess.STDOUT, cwd=out, env=env)
        manifest['pid'] = proc.pid
        mp.write_text(json.dumps(manifest, indent=2))
        try:
            proc.wait(timeout=a.timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            manifest.update(status='timeout', returncode=proc.returncode,
                            elapsed_seconds=time.monotonic()-begin)
            mp.write_text(json.dumps(manifest, indent=2))
            raise
    manifest.update(returncode=proc.returncode, elapsed_seconds=time.monotonic()-begin)
    log = (out / 'duo.out').read_text(errors='replace')
    states = out / f'AlF_refined_J{a.jmax}_unbound.states'
    try:
        if proc.returncode or 'The transition intensities are not requested (states_only option)' not in log:
            raise RuntimeError('Duo did not complete states_only output; inspect duo.out')
        rows = read_states(states)
        expected = set(range(a.jmax+1))
        if {r['J'] for r in rows} != expected:
            raise RuntimeError('Incomplete J coverage')
        zpes = set(re.findall(r'(?i)Zero point energy \(ZPE\)\s*=\s*([-\d.]+)', log))
        numeric_zpes = [float(x) for x in zpes]
        if not numeric_zpes or max(numeric_zpes)-min(numeric_zpes)>1e-6:
            raise ValueError(f'Expected one consistent energy zero; found {zpes}')
        manifest.update(status='complete', states_count=len(rows), zpe_cm=float(max(zpes, key=len)),
                        states_sha256=hashlib.sha256(states.read_bytes()).hexdigest())
    except Exception as e:
        manifest.update(status='failed', error=str(e))
        mp.write_text(json.dumps(manifest, indent=2))
        raise
    try:
        summarize(source.read_text(), rows, manifest['zpe_cm'], out, dissociation=a.dissociation)
    except Exception as e:
        manifest.update(status='analysis_failed',error=str(e))
        mp.write_text(json.dumps(manifest, indent=2))
        raise
    mp.write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
