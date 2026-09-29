"""Long-box convergence of isolated A(f) outer-well states only.

A Dirichlet wall at 3 A lies deep in the barrier; these outer levels have
negligible amplitude there. This check does not include X/A e-parity coupling
and must not be used as a count for the complete coupled molecular system.
"""
import argparse
import json
from barrier_spectrum import KINETIC, LIMIT, ZPE, effective, write_csv
import numpy as np
from scipy.linalg import eigh_tridiagonal
from coupled import fields
from pathlib import Path


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    f=fields(a.input.read_text())['POTEN',2,2]
    rows=[]
    for rmax,h in [(80.,.005),(160.,.005),(160.,.0025)]:
        r=np.arange(3.+h,rmax,h)
        off=np.full(len(r)-1,-KINETIC/h**2)
        for j in range(1,19):
            diag=2*KINETIC/h**2+effective(r,f,j)-LIMIT
            e,p=eigh_tridiagonal(diag,off,select='v',select_range=(-25.,0.),tol=1e-10,
                                 check_finite=False)
            for n,energy in enumerate(e):
                rows.append(dict(rmax=rmax,step=h,J=j,outer_index=n,binding_cm=float(-energy),
                                 mean_r=float(r@(p[:,n]**2))))
        print(rmax,h,'finished',flush=True)
    write_csv(a.output/'outer_well_convergence.csv',rows)
    selected=[x for x in rows if x['rmax']==160 and x['step']==.0025]
    (a.output/'summary.json').write_text(json.dumps(dict(f_outer_bound_count=len(selected),
        J1=[x for x in selected if x['J']==1],max_J=max(x['J'] for x in selected)),indent=2))


if __name__=='__main__':main()
