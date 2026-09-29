"""Assemble the AlF native-Duo and independent radial diagnostics."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import Counter
from barrier_spectrum import np, effective, inner_barrier, LIMIT, ZPE, write_csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bound_states import read_states
from coupled import fields, observations


def read_csv(path):
    with Path(path).open() as f:return list(csv.DictReader(f))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('work',type=Path);ap.add_argument('--model',type=Path,required=True)
    ap.add_argument('--observations',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    def native(name,jmax):return read_states(a.work/name/f'AlF_refined_J{jmax}_unbound.states')
    old=native('J160_v100',160);new=native('J160_v180',160)
    n8=native('J1_box8_v250',1);n12=native('J1_box12_v400',1)
    dvr=read_csv(a.work/'dvr601/spectrum.csv')
    checked=read_csv(a.work/'dvr801_high/spectrum.csv')
    limit=LIMIT-ZPE
    def bound(rows):return [x for x in rows if x['state']=='A1Pi' and x['flag']=='b' and x['energy']<limit]
    b=bound(new)
    def key(r):return r['J'],r['ef'],r['v']
    ref={key(x):x for x in bound(old)}
    assert set(ref)=={key(x) for x in b}
    delta=max(abs(x['energy']-ref[key(x)]['energy']) for x in b)
    obs=[x for x in observations(a.observations.read_text()) if x['state']==2 and x['J']>=0 and x['observed_cm']>=0]
    lookup={(x['J'],x['parity'],x['v']):x for x in b}
    matches=[]
    for o in obs:
        c=lookup[o['J'],o['parity'],o['v']]
        matches.append(dict(J=o['J'],parity=o['parity'],v=o['v'],observed_cm=o['observed_cm'],
                            calculated_cm=c['energy'],flag=c['flag'],mean_r=c['mean_r'],tail_probability=c['tail_probability'],
                            binding_cm=limit-c['energy']))
    vmax=[]
    for v in range(17):
        ov=[x for x in obs if x['v']==v];bv=[x for x in b if x['v']==v]
        vmax.append(dict(v=v,observed_Jmax=max((x['J'] for x in ov),default=''),
                         bound_Jmax_e=max(x['J'] for x in bv if x['ef']=='e'),
                         bound_Jmax_f=max(x['J'] for x in bv if x['ef']=='f'),
                         observed_rows=len(ov)))
    write_csv(a.output/'observed_vs_bound.csv',vmax)
    write_csv(a.output/'observed_classifications.csv',matches)
    write_csv(a.output/'A_inner_bound_levels.csv',b)
    write_csv(a.output/'A_J1_native_box8.csv',[x for x in n8 if x['state']=='A1Pi'])
    write_csv(a.output/'A_J1_native_box12.csv',[x for x in n12 if x['state']=='A1Pi'])
    radial={(int(x['J']),int(x['inner_v'])):float(x['energy']) for x in dvr
            if x['inner_v']!='' and x['below_limit']=='True' and float(x['inner_probability'])>.5}
    dvrdiff=max(abs(x['energy']-radial[x['J'],x['v']]) for x in b if x['ef']=='f')
    griddiff=max(abs(float(x['energy'])-radial[int(x['J']),int(x['inner_v'])]) for x in checked
                 if x['inner_v']!='' and x['below_limit']=='True' and float(x['inner_probability'])>.5)
    sensitivity=[]
    for name,data in [('box8_v250',n8),('box12_v400',n12)]:
        f=[x for x in data if x['state']=='A1Pi' and x['ef']=='f']
        for t in [.1,.01,.001,.0001]:
            sensitivity.append(dict(run=name,tail_threshold=t,mean_r_threshold=4,
                below_limit_localized=sum(x['energy']<limit and x['mean_r']<=4 and x['tail_probability']<=t for x in f),
                above_limit_localized=sum(x['energy']>=limit and x['mean_r']<=4 and x['tail_probability']<=t for x in f)))
    write_csv(a.output/'threshold_sensitivity.csv',sensitivity)
    fs=fields(a.model.read_text());f=fs['POTEN',2,2]
    barriers=[dict(J=j,barrier_r=inner_barrier(f,j)[0],barrier_energy=inner_barrier(f,j)[1]-ZPE) for j in [1,40,92,144,145,200,260,280]]
    outer=json.loads((a.work/'outer_check/summary.json').read_text())
    summary=dict(dissociation_cm=limit,zpe_cm=ZPE,inner_bound_count=len(b),parity_counts=dict(Counter(x['ef'] for x in b)),
                 lowest_J=1,lowest_J_inner_vmax=16,highest_inner_bound_J=144,
                 observed_vmax=max(x['v'] for x in obs),observed_Jmax=max(x['J'] for x in obs),
                 observed_raw_rows_checked=len(matches),observed_positive_weight_rows=sum(x['input_weight']>0 for x in obs),
                 all_observations_localized_below_limit=True,
                 highest_observed_energy=max(x['observed_cm'] for x in obs),
                 observed_min_binding_cm=min(x['binding_cm'] for x in matches),
                 basis_100_to_180_bound_max_delta_cm=delta,
                 independent_sinc_601_vs_native180_bound_max_delta_cm=dvrdiff,
                 independent_601_to_801_selected_bound_max_delta_cm=griddiff,
                 barriers=barriers,outer_f_channel=outer,
                 scientific_status='Extrapolation beyond observed v/J; localized above-limit states are resonance candidates only.')
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2))
    manifests={}
    for name in ['J40_native','J160_v100','J160_v180','J1_box8_v250','J1_box12_v400']:
        m=json.loads((a.work/name/'manifest.json').read_text());assert m['status']=='complete'
        manifests[name]=m
    (a.output/'provenance.json').write_text(json.dumps(dict(model=str(a.model.resolve()),
        model_sha256=hashlib.sha256(a.model.read_bytes()).hexdigest(),observations=str(a.observations.resolve()),
        observations_sha256=hashlib.sha256(a.observations.read_bytes()).hexdigest(),native_runs=manifests),indent=2))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(12,4.9),layout='constrained')
    ax=axs[0]
    js=np.arange(1,261)
    pe=np.array([inner_barrier(f,int(j))[1]-ZPE for j in js])
    ax.plot(js,pe/1000,color='#a34c17',label='A(f) effective barrier')
    for v,col in [(0,'#23724d'),(6,'#336f9f'),(16,'#7962a5')]:
        pts=[x for x in dvr if x['inner_v']==str(v) and float(x['inner_probability'])>.9]
        ax.plot([int(x['J']) for x in pts],[float(x['energy'])/1000 for x in pts],color=col,label=f'Inner sequence v = {v}')
    ax.axhline(limit/1000,color='#333',ls='--',lw=1.2,label='Dissociation limit')
    ax.scatter([x['J'] for x in obs],[x['observed_cm']/1000 for x in obs],s=5,color='black',alpha=.4,label='Observed A levels',rasterized=True)
    ax.set(xlim=(0,265),ylim=(43,80),xlabel='J',ylabel='Energy above X(v=0, J=0) / 1000 cm⁻¹',title='Binding and the rotational barrier')
    ax.legend(fontsize=8,loc='upper left')
    ax=axs[1]
    ax.barh([x['v'] for x in vmax],[x['bound_Jmax_f'] for x in vmax],height=.85,color='#c7e5d7',label='Inner-well levels below dissociation')
    ax.scatter([x['J'] for x in obs],[x['v'] for x in obs],s=5,color='black',alpha=.45,label='Observed levels',rasterized=True)
    ax.set(xlim=(0,155),ylim=(-.7,16.7),xlabel='J',ylabel='v',yticks=range(17),title='Observed coverage is well inside the bound domain')
    ax.legend(fontsize=8,loc='upper right')
    fig.savefig(a.output/'bound_levels.png',dpi=190);fig.savefig(a.output/'bound_levels.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    r=np.linspace(1.15,8,3500);V=effective(r,f,1)-ZPE
    axs[0].plot(r,V,color='#294a6c')
    axs[0].axhline(limit,ls='--',color='#444',lw=1)
    for row in sorted([x for x in n8 if x['state']=='A1Pi' and x['ef']=='f' and x['flag']=='b'],key=lambda x:x['energy']):
        axs[0].hlines(row['energy'],1.45,2.25,color='#23724d' if row['energy']<limit else '#c47719',lw=.8)
    axs[0].set(xlim=(1.2,4),ylim=(43200,62000),xlabel='r / Å',ylabel='Energy / cm⁻¹',title='J = 1: 17 inner bound + localized candidates')
    r=np.linspace(3.5,35,2000)
    axs[1].plot(r,effective(r,f,1)-LIMIT,color='#294a6c',label='Outer well, A(f), J = 1')
    for row in outer['J1']:
        axs[1].hlines(-row['binding_cm'],4,row['mean_r']+2,color='#337dab',lw=1)
    axs[1].axhline(0,color='#444',ls='--',lw=1)
    axs[1].axvline(4,color='#a34c17',ls=':',label='Mean-r cutoff = 4 Å')
    axs[1].set(xlim=(3.5,30),ylim=(-19,3),xlabel='r / Å',ylabel='Energy relative to dissociation / cm⁻¹',title='Outer bound levels fail the localization cutoff')
    axs[1].legend(fontsize=8)
    fig.savefig(a.output/'inner_and_outer_wells.png',dpi=190);fig.savefig(a.output/'inner_and_outer_wells.pdf');plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
