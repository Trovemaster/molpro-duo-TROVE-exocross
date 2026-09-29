"""Constrained A(f) surrogate refinement, always requiring native validation.

All original positive-weight A observations participate. Native e/f differences
are frozen at the supplied baseline, and are rechecked with the full X/A Duo
Hamiltonian after fitting. The output is a trial, never a production selection.
"""
import argparse
import csv
import json
import hashlib
import os
from pathlib import Path
import sys
import time
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[n]='1'
deps=Path(__file__).resolve().parent.parent/'.pecfit_stage/plot_deps'
if deps.exists():sys.path.insert(0,str(deps))
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares
from coupled import fields
from barrier import components, extrema
from barrier_spectrum import kinetic, KINETIC


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('baseline',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--height',type=float,required=True);ap.add_argument('--nfev',type=int,default=250)
    ap.add_argument('--basis',type=int,default=80)
    ap.add_argument('--loss-scale',type=float,default=.1)
    ap.add_argument('--linear',action='store_true')
    ap.add_argument('--start',type=Path)
    ap.add_argument('--flexible-crossing',action='store_true')
    ap.add_argument('--central-step',type=float,default=1e-4)
    ap.add_argument('--jacobian',choices=['central','relative'],default='central')
    ap.add_argument('--minimum-reference',action='store_true',
                    help='Constrain height above the actual adiabatic minimum; otherwise use V(RE_EMO).')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    text=(a.baseline/'input.inp').read_text();f=fields(text)['POTEN',2,2];p=np.array(f.values)
    if f.kind!='COUPLED-PEC' or len(p)!=41 or tuple(p[4:8])!=(6.,6.,6.,6.):
        raise ValueError('This experimental AlF surrogate requires the 15/12/13, PL=PR=6 layout')
    if a.loss_scale<=0 or a.central_step<=0 or a.basis<8:
        raise ValueError('Positive loss/derivative scales and at least eight basis functions required')
    (a.output/'manifest.json').write_text(json.dumps(dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
        baseline_sha256=hashlib.sha256((a.baseline/'input.inp').read_bytes()).hexdigest(),
        residuals_sha256=hashlib.sha256((a.baseline/'residuals.csv').read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
    (a.output/'fit_low_barrier_snapshot.py').write_bytes(Path(__file__).read_bytes())
    rows=list(csv.DictReader((a.baseline/'residuals.csv').open()))
    obs=[x for x in rows if x['state']=='2' and float(x['input_weight'])>0]
    js=sorted({int(float(x['J'])) for x in obs});jmap={j:i for i,j in enumerate(js)}
    ji=np.array([jmap[int(float(x['J']))] for x in obs]);vi=np.array([int(x['v']) for x in obs])
    weights=np.sqrt([float(x['input_weight']) for x in obs])
    r=np.linspace(.7,8,601);T=kinetic(r)
    def V(rr,q):
        u,v,w=components(np.asarray(rr),q);return (u+v-np.hypot(u-v,2*w))*.5
    original=V(r,p);H=T.copy();H.flat[::len(r)+1]+=original
    ev,U=eigh(H,subset_by_index=(0,a.basis-1),check_finite=False)
    Br=U.T@((KINETIC/r**2)[:,None]*U)
    def eigenvalues(q):
        h=np.diag(ev)+U.T@((V(r,q)-original)[:,None]*U)
        return np.array([eigh(h+(j*(j+1)-2)*Br,subset_by_index=(0,6),eigvals_only=True,
                              check_finite=False,driver='evr') for j in js])
    old=eigenvalues(p)
    target=old[ji,vi]+np.array([float(x['residual_cm']) for x in obs])
    def unpack(x):
        q=p.copy();q[0]=1000*x[0];q[1]=x[1];q[8:15]=x[2:9]
        if a.flexible_crossing:
            q[16]=np.exp(x[9]);q[27]=1000*x[10];q[28]=x[11];q[35]=x[12]
        return q
    rr=np.linspace(.7,12,241)
    ps=np.array(fields(a.start.read_text())['POTEN',2,2].values) if a.start else p
    x0=np.r_[ps[0]/1000,ps[1],ps[8:15]]
    lo=np.r_[42,1.60,.3,np.full(6,-30)]
    hi=np.r_[47,1.70,5.,np.full(6,30)]
    if a.flexible_crossing:
        x0=np.r_[x0,np.log(ps[16]),ps[27]/1000,ps[28],ps[35]]
        lo=np.r_[lo,np.log(1e4),.01,1.9,.1]
        hi=np.r_[hi,np.log(1e9),10.,3.,3.]
    ex=extrema(r,V(r,ps));rb=next(x[0] for x in ex if x[2]=='max')
    x0=np.r_[x0,rb];lo=np.r_[lo,1.97];hi=np.r_[hi,2.35]
    calls=0;begin=time.monotonic()
    def fun(x):
        nonlocal calls
        calls+=1;q=unpack(x);delta=(eigenvalues(q)[ji,vi]-target)*weights
        # Cauchy data loss with a 0.1 cm^-1 scale, quadratic shape constraints.
        scale=a.loss_scale;rob=delta if a.linear else np.sign(delta)*scale*np.sqrt(np.log1p((delta/scale)**2))
        re_=q[1];y=(rr**6-re_**6)/(rr**6+re_**6);dy=12*rr**5*re_**6/(rr**6+re_**6)**2
        beta=np.polynomial.polynomial.polyval(y,q[8:15]);db=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(q[8:15]))*dy
        rb=x[-1];vv=V([rb-.001,rb,rb+.001],q);minimum=V(q[1],q)
        if a.minimum_reference:
            from scipy.optimize import minimize_scalar
            minimum=minimize_scalar(lambda rr:float(V(rr,q)),bounds=(1.55,1.75),method='bounded',
                                    options={'xatol':1e-10}).fun
        shape=np.r_[(vv[1]-minimum-a.height)/.5,(vv[2]-vv[0])/.002/30,
                    max(0,(vv[0]+vv[2]-2*vv[1])/.001**2+500)/100,
                    np.maximum(.05-beta,0)*1000,np.maximum(.05-beta-(rr-re_)*db,0)*1000]
        if calls%150==0:print(json.dumps(dict(calls=calls,elapsed=time.monotonic()-begin,
            rms=float(np.sqrt(np.mean(delta**2))),median=float(np.median(abs(delta))),
            barrier=float(vv[1]),radius=float(rb),shape_norm=float(np.linalg.norm(shape)))),flush=True)
        return np.r_[rob,shape,1e-4*(x[2:9]-x0[2:9])]
    def jac(x):
        # Absolute, symmetric steps avoid round-off in the curvature constraint
        # and in small high-order EMO coefficients near zero.
        columns=[]
        for i in range(len(x)):
            step=a.central_step;xp=x.copy();xm=x.copy();xp[i]+=step;xm[i]-=step
            columns.append((fun(xp)-fun(xm))/(2*step))
        return np.array(columns).T
    result=least_squares(fun,x0,bounds=(lo,hi),x_scale='jac',max_nfev=a.nfev,
                         jac=jac if a.jacobian=='central' else '2-point',diff_step=2e-6,
                         ftol=1e-9,xtol=1e-10,gtol=1e-8)
    q=unpack(result.x)
    for span,value in reversed(list(zip(f.spans,q))):
        start,stop=span;text=text[:start]+f'{value:.14E}'+text[stop:]
    (a.output/'trial.inp').write_text(text)
    delta=eigenvalues(q)[ji,vi]-target
    rg=np.unique(np.r_[np.linspace(.7,8,10001),np.geomspace(8,1000,1000)])
    info=dict(success=bool(result.success),message=result.message,nfev=result.nfev,calls=calls,
        elapsed_seconds=time.monotonic()-begin,target_height=a.height,
        linear_loss=a.linear,loss_scale=a.loss_scale,flexible_crossing=a.flexible_crossing,
        height_reference='adiabatic minimum' if a.minimum_reference else 'V(RE_EMO)',
        jacobian=a.jacobian,central_step=a.central_step,
        surrogate_rms_cm=float(np.sqrt(np.mean(delta**2))),surrogate_median_cm=float(np.median(abs(delta))),
        extrema=extrema(rg,V(rg,q)),basis=a.basis,
        warning='Surrogate result: full native Duo validation and basis checks required.')
    (a.output/'result.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))


if __name__=='__main__':main()
