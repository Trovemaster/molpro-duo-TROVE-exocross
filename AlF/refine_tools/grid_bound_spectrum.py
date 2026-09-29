"""Resolve inner/outer character using the actual PEC grid printed by native Duo.

This applies to the AlF A(f) channel without additional state-2 operators.
It never infers vibrational character from the native localization flag alone.
"""
import argparse
import csv
import json
import os
from pathlib import Path
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='2'
import numpy as np
from scipy.linalg import eigh
from barrier_spectrum import kinetic,KINETIC,write_csv
from bound_states import read_states


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run',type=Path);ap.add_argument('--js',nargs='+',type=int,default=[1])
    ap.add_argument('--dissociation',type=float,default=55564.02471)
    a=ap.parse_args();manifest=json.loads((a.run/'manifest.json').read_text())
    if manifest.get('status')!='complete':raise ValueError('Native run must be complete')
    zpe=manifest['zpe_cm'];data=np.loadtxt(a.run/'Potential_functions.dat');r=data[:,0];v=data[:,2]
    if not np.allclose(np.diff(r),np.diff(r).mean(),rtol=1e-7,atol=1e-9):raise ValueError('Uniform DVR required')
    # Printed radii are rounded to ten decimals; reconstruct the exact uniform grid.
    r=np.linspace(r[0],r[-1],len(r))
    T=kinetic(r);states=read_states(next(a.run.glob('*.states')));summary=[];levels=[]
    for j in a.js:
        veff=v+KINETIC*(j*(j+1)-2)/r**2
        ip=np.flatnonzero((veff[1:-1]>veff[:-2])&(veff[1:-1]>veff[2:]))+1
        ip=[i for i in ip if 2.0<r[i]<3.5]
        if len(ip)!=1:raise ValueError(f'Expected one inner barrier: J={j}, indices={ip}')
        rb=r[ip[0]];peak=veff[ip[0]]
        H=T.copy();H.flat[::len(r)+1]+=veff
        e,u=eigh(H,subset_by_value=(42000,max(a.dissociation,peak)+200),check_finite=False,driver='evr')
        pin=np.sum(u[r<rb]**2,axis=0);means=r@(u*u)
        native={x['v']:x for x in states if x['state']=='A1Pi' and x['J']==j and x['ef']=='f'}
        errs=[];inner=[]
        for k,energy in enumerate(e):
            n=native.get(k)
            # Comparison by ordered f-parity eigenvalue, including outer/box levels.
            error=None if n is None else float(energy-zpe-n['energy'])
            if pin[k]>.5 and energy<a.dissociation:inner.append(k)
            if k<7 and error is not None:errs.append(abs(error))
            levels.append(dict(J=j,index=k,energy=float(energy-zpe),inner_probability=float(pin[k]),mean_r=float(means[k]),
                below_limit=bool(energy<a.dissociation),below_barrier=bool(energy<peak),native_flag=n['flag'] if n else '',
                native_error_cm=error))
        summary.append(dict(J=j,inner_bound_count=len(inner),inner_vmax=len(inner)-1,
            total_below_limit=int(np.count_nonzero(e<a.dissociation)),
            outer_below_limit=int(np.count_nonzero((e<a.dissociation)&(pin<=.5))),
            inner_below_barrier=int(np.count_nonzero((e<peak)&(pin>.5))),barrier_r=float(rb),barrier_cm=float(peak),
            minimum_inner_probability=min((float(pin[k]) for k in inner),default=None),
            observed_v0_v6_max_native_error_cm=max(errs,default=None)))
    write_csv(a.run/'inner_outer_levels.csv',levels)
    (a.run/'inner_outer_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
