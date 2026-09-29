"""Extract a named PEC from Qin's independently paired radial columns.

Each state has its own radius column. Missing trailing pairs are not zeroes,
and the X-state radii must never be reused for the A-state energies.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def extract_pair(text, state='A1Pi'):
    lines = text.splitlines()
    header = lines[0].split()
    if len(header) % 2 or any(x != 'R/Angstrom' for x in header[::2]):
        raise ValueError('Expected alternating R/Angstrom and state-name columns')
    if header[1::2].count(state) != 1:
        raise ValueError(f'Expected exactly one {state} energy column')
    col = header.index(state)-1
    rows = []
    for lineno, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        # Preserve empty cells so another state's values cannot slide left.
        cells = line.split('\t') if '\t' in line else line.split()
        cells = [s.strip() for s in cells]
        if len(cells) > len(header) and any(cells[len(header):]):
            raise ValueError(f'Unexpected extra columns on line {lineno}')
        cells += [''] * max(0, len(header)-len(cells))
        x, y = cells[col:col+2]
        if not x and not y:
            continue
        if not x or not y:
            raise ValueError(f'Incomplete {state} pair on line {lineno}')
        float(x); float(y)
        rows.append((x, y, lineno))
    if not rows:
        raise ValueError('No reference points found')
    return rows, [col+1, col+2]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source', type=Path)
    ap.add_argument('model', type=Path)
    ap.add_argument('output', type=Path, help='New package directory')
    ap.add_argument('--fit-factor', type=float, default=100)
    ap.add_argument('--energy-unit', required=True, choices=['cm-1', 'hartree'])
    a = ap.parse_args()
    raw = a.source.read_bytes()
    rows, columns = extract_pair(raw.decode('utf-8-sig'))
    a.output.mkdir(parents=True, exist_ok=False)
    (a.output/'Potential_energy_original.txt').write_bytes(raw)
    data = a.output/'A_original_angstrom_cm-1.dat'
    if a.energy_unit != 'cm-1':
        data = a.output/'A_original_angstrom_hartree.dat'
    data.write_text('# Qin et al. 2022, supplied source A1Pi pair, original values\n'
                    f'# r_angstrom V_{a.energy_unit}\n'+
                    ''.join(f'{x} {y}\n' for x, y, _ in rows))
    metadata = dict(source_filename=a.source.name, source_sha256=hashlib.sha256(raw).hexdigest(),
                    source_columns_one_based=columns, state='A1Pi', point_count=len(rows),
                    source_line_numbers=[n for _, _, n in rows],
                    radius_unit='angstrom', energy_unit=a.energy_unit,
                    energy_unit_note='Explicit caller choice; the source labels radius units only.',
                    citation='Qin, Bai and Liu 2022, MNRAS 510, 3011, DOI 10.1093/mnras/stab3598')
    (a.output/'extraction.json').write_text(json.dumps(metadata, indent=2)+'\n')
    command = [sys.executable, str(Path(__file__).with_name('add_abinitio_reference.py')),
               str(a.model), str(data), str(a.output/'inputs'), '--distance-unit', 'angstrom',
               '--energy-unit', a.energy_unit, '--fit-factor', str(a.fit_factor),
               '--citation', metadata['citation']]
    subprocess.run(command, check=True)
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
