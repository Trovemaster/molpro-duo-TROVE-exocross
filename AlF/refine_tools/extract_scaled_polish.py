"""Extract the last complete native proposal, never its preceding residual table."""
from pathlib import Path
import argparse
import json
from coupled import checkpoints,transfer
from barrier import audit
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('input',type=Path);ap.add_argument('log',type=Path);ap.add_argument('output',type=Path)
a=ap.parse_args()
source=a.input.read_text()
cp=checkpoints(a.log.read_text())
if not cp:raise ValueError('No complete native proposal')
text=transfer(source,cp[-1])
shape=audit(text,dict(barrier_r_min=2.2,barrier_r_max=3.1,barrier_above_limit_min=100,barrier_above_limit_max=8000,outer_well_depth_max=100))
if not shape['ok']:raise ValueError(shape)
a.output.write_text(text)
a.output.with_suffix('.shape.json').write_text(json.dumps(shape,indent=2)+'\n')
print(json.dumps(dict(checkpoints=len(cp),shape=shape),indent=2))
