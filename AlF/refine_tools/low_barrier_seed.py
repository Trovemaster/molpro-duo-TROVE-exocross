"""Test lower A-state barrier constraints while retaining the Al+F asymptote."""
import argparse
import json
import os
from pathlib import Path
import re
import sys
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
deps=Path(__file__).resolve().parent.parent/'.pecfit_stage/plot_deps'
if deps.exists():sys.path.insert(0,str(deps))
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from coupled import fields
from barrier import components, potential, extrema


def build(source,out,height,reference='asymptote'):
    text=source.read_text();fs=fields(text);f=fs['POTEN',2,2]
    p=np.array(f.values);limit=p[15]
    r=np.linspace(1.30,1.96,201);target=potential(r,f)
    minimum=minimize_scalar(lambda x:float(potential(np.array(x),f)),bounds=(1.5,1.8),method='bounded').fun
    goal=(limit if reference=='asymptote' else minimum)+height
    def unpack(x):
        q=p.copy();q[0]=x[0]*1000;q[1]=x[1];q[8:15]=x[2:9];q[16]=np.exp(x[9])
        return q
    def V(r,q):
        a,b,w=components(np.array(r),q)
        return (a+b-np.hypot(a-b,2*w))*.5
    x0=np.r_[p[0]/1000,p[1],p[8:15],np.log(p[16]),2.4]
    lo=np.r_[40,1.60,.2,np.full(6,-30),np.log(1e4),1.98]
    hi=np.r_[50,1.70,5.,np.full(6,30),np.log(1e9),2.9]
    def residual(x):
        q=unpack(x);rb=x[10];v=V([rb-.002,rb,rb+.002],q)
        rr=np.linspace(.7,12,201);re_=q[1];y=(rr**6-re_**6)/(rr**6+re_**6)
        dy=12*rr**5*re_**6/(rr**6+re_**6)**2
        beta=np.polynomial.polynomial.polyval(y,q[8:15])
        db=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(q[8:15]))*dy
        exponent_slope=beta+(rr-re_)*db
        return np.r_[(V(r,q)-target)/(1+(target-minimum)/2500),
                     (v[1]-goal)/.2,(v[2]-v[0])/.004/10,
                     np.maximum((v[0]+v[2]-2*v[1])/.002**2,0)/100,
                     np.maximum(.05-beta,0)*10000,np.maximum(.05-exponent_slope,0)*10000,
                     .01*(x[2:9]-x0[2:9])]
    best=None
    for radius in [2.35,2.10,2.6]:
        start=x0.copy();start[10]=radius
        result=least_squares(residual,start,bounds=(lo,hi),x_scale='jac',max_nfev=1800,
                             ftol=1e-11,xtol=1e-11,gtol=1e-9)
        if best is None or sum(result.fun**2)<sum(best.fun**2):best=result
        if np.max(np.abs(result.fun[:201]))<1 and abs(result.fun[201])<1:break
    q=unpack(best.x)
    # Keep all original physical parameters except state-2 EMO coefficients and
    # the repulsive amplitude. Freeze X and its BOB and fit only A plus Lx.
    edits=[]
    for key,field in fs.items():
        for i,(span,value) in enumerate(zip(field.spans,field.values)):
            start,stop=span;end=text.find('\n',stop);suffix=text[stop:end]
            suffix=re.sub(r'\bfit\b','',suffix,flags=re.I)
            fitting=(key==('POTEN',2,2) and i in [0,1,*range(8,15)]) or (key==('L+',2,1) and i==4)
            val=q[i] if key==('POTEN',2,2) else value
            edits.append((start,end,f'{val:.14E}'+(' fit' if fitting else '')+suffix))
    for start,end,value in sorted(edits,reverse=True):text=text[:start]+value+text[end:]
    text=re.sub(r'(?im)^\s*npoints\s+.*$', ' npoints 601',text,count=1)
    text=re.sub(r'(?im)^\s*range\s+.*$', ' range 0.7 8',text,count=1)
    text=re.sub(r'(?im)^\s*vmax\s+.*$', ' vmax 100 100',text,count=1)
    out.mkdir(parents=True,exist_ok=False);(out/'seed.inp').write_text(text)
    rr=np.unique(np.r_[np.linspace(.7,8,10001),np.geomspace(8,1000,1000)])
    v=V(rr,q);ex=extrema(rr,v)
    info=dict(target_reference=reference,target_height_cm=height,target_absolute_cm=goal,
              attained_target_value=float(V(best.x[10],q)),target_radius=float(best.x[10]),
              seed_curve_rms_cm=float(np.sqrt(np.mean((V(r,q)-target)**2))),
              seed_curve_max_cm=float(np.max(np.abs(V(r,q)-target))),
              optimizer_success=bool(best.success),nfev=best.nfev,extrema=ex,
              finite_on_validation_grid=bool(np.all(np.isfinite(v))),
              asymptote_absolute_cm=limit,minimum_original_absolute_cm=float(minimum),
              note='This is a curve-initialization trial, not an experimental fit.')
    (out/'seed.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2),flush=True)
    return info


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--height',type=float,required=True);ap.add_argument('--reference',choices=['asymptote','minimum'],default='asymptote')
    a=ap.parse_args();build(a.input,a.output,a.height,a.reference)
