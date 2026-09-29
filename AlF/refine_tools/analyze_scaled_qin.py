"""Report scaled-Qin grid and constrained-fit experiments using fresh native energies."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from coupled import fields,observations,save_assessment
from barrier import potential
from bound_states import read_states
from add_abinitio_reference import model_minimum,ps1997_model_weights


def stats(rows,state):
    a=[x for x in rows if int(x['state'])==state and float(x['input_weight'])>0]
    d=np.array([float(x['residual_cm']) for x in a]);w=np.array([float(x['input_weight']) for x in a])
    return dict(n=len(a),weighted_rms=float(np.sqrt(np.sum(w*d*d)/w.sum())),median_abs=float(np.median(abs(d))),max_abs=float(max(abs(d))))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('root',type=Path)
    a=ap.parse_args();root=a.root.resolve();work=root.parent
    baseline=work/'qin_2022/trial_A_only/baseline'
    original=(baseline/'input.inp').read_text();te,_=model_minimum(fields(original)['POTEN',2,2]);limit=55564.02471
    reference=np.loadtxt(root/'inputs/A_scaled_nominal_zero.dat');wp=ps1997_model_weights(original,fields(original)['POTEN',2,2],reference[:,0],.001,18000)
    rowsets={};models={}
    runs=dict(selected=baseline,weak_full=root/'final_soft',strong_full=root/'final_strong',outer=root/'final_outer')
    for name,run in runs.items():
        if name!='selected':save_assessment(run,98)
        rows=list(csv.DictReader((run/'residuals.csv').open()));rowsets[name]=rows
        f=fields((run/'input.inp').read_text())['POTEN',2,2]
        minimum,re_=model_minimum(f)
        peak=minimize_scalar(lambda r:-float(potential(np.array([r]),f)[0]),bounds=(2.2,3.1),method='bounded')
        dv=potential(reference[:,0],f)-reference[:,1]
        mask=reference[:,0]>=2.2
        models[name]=dict(A=stats(rows,2),X=stats(rows,1),Te=minimum,Re=re_,barrier_r=float(peak.x),
            barrier_above_limit=float(-peak.fun-limit),barrier_above_minimum=float(-peak.fun-minimum),
            PEC_rms_all=float(np.sqrt(np.sum(wp*dv**2)/wp.sum())),
            PEC_rms_outer=float(np.sqrt(np.sum(wp[mask]*dv[mask]**2)/wp[mask].sum())),
            input_sha256=hashlib.sha256((run/'input.inp').read_bytes()).hexdigest())
        assert json.loads((run/'assessment.json').read_text())['assignment_mismatches']==0
    gridrun=root/'grid_nominal_J98';states=read_states(next(gridrun.glob('*.states')))
    lookup={(x['J'],x['parity'],1 if x['state']=='X1Sigma+' else 2,x['v']):x for x in states}
    gridrows=[]
    for o in observations(original):
        c=lookup[(int(o['J']),o['parity'],o['state'],o['v'])]
        gridrows.append(dict(o,calculated_cm=c['energy'],residual_cm=o['observed_cm']-c['energy'],flag=c['flag']))
    rowsets['grid']=gridrows
    models['grid']=dict(A=stats(gridrows,2),X=stats(gridrows,1),energy_precision='six-decimal native states')
    with (root/'grid_residuals.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(gridrows[0]));writer.writeheader();writer.writerows(gridrows)
    grid_spectra={name:json.loads((root/name/'inner_outer_summary.json').read_text()) for name in
        ['grid_nominal_J98','grid_nominal_box9','grid_endpoint_box9','bound_soft_J1','bound_soft_box9','bound_outer_final_J1']}
    def inner_energies(name):
        ss=read_states(next((root/name).glob('*.states')))
        return np.array([x['energy'] for x in ss if x['state']=='A1Pi' and x['J']==1 and x['ef']=='f'])[:17]
    checks=dict(grid_box8_to_box9_inner_max_change_cm=float(max(abs(inner_energies('grid_nominal_J98')-inner_energies('grid_nominal_box9')))),
                soft_box8_to_box9_inner_max_change_cm=float(max(abs(inner_energies('bound_soft_J1')-inner_energies('bound_soft_box9')))))
    result=dict(scaling=json.loads((root/'inputs/scaling.json').read_text()),models=models,spectra=grid_spectra,convergence=checks,
        selected_model_replaced=False,stats_convention='Same 1834 merged input records, 1822 positive weights: 731 X and 1091 A; 12 pre-existing zero weights retained.')
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    colors=dict(selected='#7b8189',weak_full='#087f70',outer='#1b65b9',grid='#8c45a8')
    labels=dict(selected='Current selected',weak_full='Weak full-curve constraint',outer='Outer-region constraint',grid='Scaled Qin grid')
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(12.5,5),layout='constrained')
    rr=np.linspace(1.3,9,3000)
    for name in ['selected','weak_full','outer']:
        f=fields((runs[name]/'input.inp').read_text())['POTEN',2,2];v=potential(rr,f)
        ax[0].plot(rr,v-te,color=colors[name],label=labels[name],lw=1.7)
        ax[1].plot(rr,v-limit,color=colors[name],lw=1.7)
    dense=np.loadtxt(gridrun/'Potential_functions.dat')
    ax[0].plot(dense[:,0],dense[:,2]-te,color=colors['grid'],label=labels['grid'],lw=1.8)
    ax[1].plot(reference[:,0],reference[:,1]-limit,'o-',ms=3,color=colors['grid'],lw=1)
    ax[0].axhline(limit-te,color='black',ls=':',lw=1);ax[1].axhline(0,color='black',ls=':',lw=1)
    ax[0].set(xlim=(1.35,4.5),ylim=(-300,18500),xlabel='r / Å',ylabel='V − adopted A minimum / cm⁻¹',title='A-state well and barrier')
    ax[1].set(xlim=(2.2,9),ylim=(-160,1700),xlabel='r / Å',ylabel='V − dissociation limit / cm⁻¹',title='Barrier and shallow outer well')
    ax[0].legend(loc='upper right',fontsize=9)
    fig.savefig(root/'PEC_comparison.png',dpi=190);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12.5,4.4),layout='constrained')
    for offset,name in enumerate(['selected','weak_full','outer']):
        med=[];p90=[]
        for v in range(7):
            d=np.array([abs(float(x['residual_cm'])) for x in rowsets[name] if int(x['state'])==2 and int(x['v'])==v and float(x['input_weight'])>0])
            med.append(np.median(d));p90.append(np.quantile(d,.9))
        ax[0].plot(range(7),med,'o-',label=labels[name],color=colors[name])
        ax[1].plot(range(7),p90,'o-',color=colors[name])
    for aa,title in zip(ax,['Median absolute residual','90th percentile absolute residual']):
        aa.set(xlabel='A-state v',ylabel='|observed − calculated| / cm⁻¹',title=title,xticks=range(7))
        aa.grid(alpha=.15)
    ax[0].legend(fontsize=9);fig.savefig(root/'residual_comparison.png',dpi=190);plt.close(fig)
    print(json.dumps(dict(models=models,convergence=checks),indent=2))


if __name__=='__main__':main()
