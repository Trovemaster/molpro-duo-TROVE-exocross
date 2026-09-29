"""Joint spectroscopic/ab-initio A-state surrogate with full native validation required.

Minimizes sum(w_exp*dE**2) + strength*sum(w_PEC*dV**2)/sum(w_PEC),
plus smooth physical shape penalties. The electronic mixing correction is
frozen only inside the optimizer. Every accepted result needs a new native
Duo calculation. X, Lx, BOB, C5/C6, asymptote and both PL=PR=6 are fixed.
"""
import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path
import time
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares
from coupled import fields
from barrier import components, audit
from barrier_spectrum import kinetic,KINETIC
from add_abinitio_reference import ps1997_model_weights,block_text


def remove_A_reference(text):
    """Reference experiments must not inherit a different A reference from baseline."""
    pattern=r'(?ims)(?:^[ \t]*\([^\n]*\)[ \t]*\n|^[ \t]*\n)*^[ \t]*abinitio[ \t]+poten[ \t]+2\b.*?^[ \t]*end[ \t]*$'
    return re.sub(pattern,'\n',text)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('baseline',type=Path);ap.add_argument('reference',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--start',type=Path);ap.add_argument('--strength',type=float,default=.001)
    ap.add_argument('--nfev',type=int,default=180);ap.add_argument('--basis',type=int,default=70)
    ap.add_argument('--loss-scale',type=float,help='Cauchy scale for experimental residuals only; omit for ordinary least squares')
    ap.add_argument('--reference-rmin',type=float,default=0.,help='Use original reference points at or above this radius')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    text=remove_A_reference((a.baseline/'input.inp').read_text());f=fields(text)['POTEN',2,2];p=np.array(f.values)
    if len(p)!=41 or tuple(p[4:8])!=(6.,6.,6.,6.):raise ValueError('Expected AlF 15/12/13 coupled PEC')
    ps=np.array(fields(a.start.read_text())['POTEN',2,2].values) if a.start else p.copy()
    reference=np.loadtxt(a.reference);rref,vref=reference.T
    wp=ps1997_model_weights(text,f,rref,.001,18000.)
    wp=np.where(rref>=a.reference_rmin,wp,0.)
    if a.strength<=0 or wp.sum()<=0 or (a.loss_scale is not None and a.loss_scale<=0):raise ValueError('Positive weights and loss scale required')
    obs=[x for x in csv.DictReader((a.baseline/'residuals.csv').open()) if x['state']=='2' and float(x['input_weight'])>0]
    js=sorted({int(float(x['J'])) for x in obs});jm={j:i for i,j in enumerate(js)}
    ji=np.array([jm[int(float(x['J']))] for x in obs]);vi=np.array([int(x['v']) for x in obs])
    we=np.array([float(x['input_weight']) for x in obs]);sw=np.sqrt(we)
    def V(rr,q):
        u,v,w=components(np.asarray(rr),q);return .5*(u+v-np.hypot(u-v,2*w))
    r=np.linspace(.7,8,601);T=kinetic(r);original=V(r,p)
    H=T.copy();H.flat[::len(r)+1]+=original
    ev,U=eigh(H,subset_by_index=(0,a.basis-1),check_finite=False)
    B=U.T@((KINETIC/r**2)[:,None]*U)
    def spectrum(q):
        h=np.diag(ev)+U.T@((V(r,q)-original)[:,None]*U)
        return np.array([eigh(h+(j*(j+1)-2)*B,subset_by_index=(0,6),eigvals_only=True,check_finite=False) for j in js])
    target=spectrum(p)[ji,vi]+np.array([float(x['residual_cm']) for x in obs])
    def pack(q):return np.r_[q[0]/1000,q[1],q[8:15],np.log(q[16]),q[27]/1000,q[28],q[35]]
    def unpack(x):
        q=p.copy();q[0]=1000*x[0];q[1]=x[1];q[8:15]=x[2:9]
        q[16]=np.exp(x[9]);q[27]=1000*x[10];q[28]=x[11];q[35]=x[12]
        return q
    lo=np.r_[42,1.60,.3,np.full(6,-20.),np.log(1e5),.01,2.,.1]
    hi=np.r_[46,1.72,5.,np.full(6,20.),np.log(1e8),12.,3.2,4.]
    x0=pack(ps);rr=np.linspace(.7,12,401);outer=np.linspace(3.2,12,181)
    begin=time.monotonic();calls=0
    def fun(x):
        nonlocal calls
        calls+=1;q=unpack(x);de=(spectrum(q)[ji,vi]-target)*sw
        dv=(V(rref,q)-vref)*np.sqrt(a.strength*wp/wp.sum())
        y=(rr**6-q[1]**6)/(rr**6+q[1]**6);dy=12*rr**5*q[1]**6/(rr**6+q[1]**6)**2
        beta=np.polynomial.polynomial.polyval(y,q[8:15]);db=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(q[8:15]))*dy
        shape=np.r_[np.maximum(.05-beta,0)*1000,np.maximum(.05-beta-(rr-q[1])*db,0)*1000,
                    np.maximum(p[15]-100-V(outer,q),0),
                    np.maximum(V(np.array([2.1,2.95]),q)-V(np.array([2.5]),q),0)]
        if calls%400==0:print(json.dumps(dict(calls=calls,seconds=time.monotonic()-begin,
            exp_rms=float(np.sqrt(np.sum(de**2)/we.sum())),PEC_rms=float(np.sqrt(np.sum(dv**2)/a.strength)),shape_norm=float(np.linalg.norm(shape)))),flush=True)
        data_loss=de if a.loss_scale is None else np.sign(de)*a.loss_scale*np.sqrt(np.log1p((de/a.loss_scale)**2))
        return np.r_[data_loss,dv,shape,1e-5*(x[2:9]-x0[2:9])]
    def jac(x):
        cols=[]
        for i in range(len(x)):
            step=1e-5;xp=x.copy();xm=x.copy();xp[i]+=step;xm[i]-=step
            cols.append((fun(xp)-fun(xm))/(2*step))
        return np.array(cols).T
    result=least_squares(fun,x0,jac=jac,bounds=(lo,hi),x_scale='jac',max_nfev=a.nfev,ftol=1e-8,xtol=1e-9,gtol=1e-7)
    q=unpack(result.x)
    for (start,stop),value in reversed(list(zip(f.spans,q))):text=text[:start]+f'{value:.14E}'+text[stop:]
    (a.output/'trial.inp').write_text(text)
    at=re.search(r'(?im)^\s*FITTING\s*$',text).start()
    block=block_text(2,'A1Pi',1,1,rref,vref,a.strength,.001,18000,wp)
    (a.output/'trial_with_reference.inp').write_text(text[:at]+'\n'+block+'\n'+text[at:])
    shape=audit(text,dict(barrier_r_min=2.2,barrier_r_max=3.1,barrier_above_limit_min=100,barrier_above_limit_max=8000,outer_well_depth_max=100))
    de=spectrum(q)[ji,vi]-target;dv=V(rref,q)-vref
    info=dict(success=bool(result.success),message=result.message,nfev=result.nfev,calls=calls,
        seconds=time.monotonic()-begin,strength=a.strength,basis=a.basis,loss_scale=a.loss_scale,reference_rmin=a.reference_rmin,
        surrogate_A_weighted_rms=float(np.sqrt(np.sum(we*de**2)/we.sum())),
        surrogate_A_median=float(np.median(abs(de))),PEC_weighted_rms=float(np.sqrt(np.sum(wp*dv**2)/wp.sum())),
        shape=shape,baseline_sha256=hashlib.sha256((a.baseline/'input.inp').read_bytes()).hexdigest(),
        source_reference_sha256=hashlib.sha256(a.reference.read_bytes()).hexdigest(),
        warning='Joint surrogate candidate; native validation required. Strength is defined by the objective in this script; do not infer native weight normalization.')
    (a.output/'result.json').write_text(json.dumps(info,indent=2)+'\n')
    (a.output/'code_snapshot.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(info,indent=2),flush=True)


if __name__=='__main__':main()
