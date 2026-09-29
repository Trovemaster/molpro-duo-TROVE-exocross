"""Uniformly scale Qin A about its sampled minimum; retain both energy-zero conventions.

No radial transformation or additional reference points are used. The primary
curve assumes the source separated-atom zero is zero cm-1. A sensitivity curve
instead treats the finite 9-A endpoint as the source asymptote. The absolute
target asymptote is the adopted model value, not a new measurement.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from coupled import fields, uncomment
from add_abinitio_reference import model_minimum, block_text, ps1997_model_weights


def scale(points, te, limit, source_limit=0.):
    points=np.asarray(points,float)
    if points.ndim!=2 or points.shape[1]!=2 or len(points)<3 or not np.isfinite(points).all():
        raise ValueError('Expected finite two-column PEC data')
    if np.any(points[:,0]<=0) or np.any(np.diff(points[:,0])<=0):
        raise ValueError('Radii must be positive and strictly increasing')
    minimum=float(points[:,1].min())
    if not np.isfinite([te,limit,source_limit]).all() or limit<=te or source_limit<=minimum:
        raise ValueError('Both well depths must be positive and finite')
    factor=(limit-te)/(source_limit-minimum)
    return np.c_[points[:,0],te+factor*(points[:,1]-minimum)],factor


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model',type=Path);ap.add_argument('data',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--dissociation',type=float,default=55564.02471)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    text=a.model.read_text();f=fields(text)['POTEN',2,2]
    te,re_=model_minimum(f);raw=np.loadtxt(a.data)
    match=re.search(r'(?ims)^poten\s+2\b.*?^end\s*$',uncomment(text))
    if not match:raise ValueError('Missing state 2')
    weights=ps1997_model_weights(text,f,raw[:,0],.001,18000.)
    output=dict(model_sha256=hashlib.sha256(a.model.read_bytes()).hexdigest(),
                data_sha256=hashlib.sha256(a.data.read_bytes()).hexdigest(),
                target_te=te,target_re=re_,target_asymptote=a.dissociation,
                target_A_depth=a.dissociation-te,reference_minimum=float(raw[:,1].min()),
                minimum_convention='sampled minimum; interpolation may lie slightly lower',curves={})
    for name,source_limit in [('nominal_zero',0.),('endpoint',float(raw[-1,1]))]:
        data,s=scale(raw,te,a.dissociation,source_limit)
        np.savetxt(a.output/f'A_scaled_{name}.dat',data,fmt=['%.12f','%.10f'],header='r_A V_cm-1; original radii, uniform energy scale')
        grid='\n'.join(['poten 2','name "A1Pi"','lambda 1','symmetry +','mult 1','type grid','units angstrom cm-1','values']+
                       [f'{r:.12f} {e:.10f}' for r,e in data]+['end'])
        (a.output/f'AlF_grid_{name}.inp').write_text(text[:match.start()]+grid+text[match.end():])
        block=block_text(2,'A1Pi',1,1,data[:,0],data[:,1],100,.001,18000,weights)
        (a.output/f'A_reference_{name}.inp').write_text(block)
        insert=re.search(r'(?im)^\s*FITTING\s*$',uncomment(text)).start()
        (a.output/f'AlF_reference_{name}.inp').write_text(text[:insert]+'\n'+block+'\n'+text[insert:])
        output['curves'][name]=dict(source_asymptote=source_limit,source_depth=source_limit-float(raw[:,1].min()),
            scale=s,sampled_barrier_cm=float(data[(raw[:,0]>=2)&(raw[:,0]<=4),1].max()),
            endpoint_cm=float(data[-1,1]))
    (a.output/'scaling.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
