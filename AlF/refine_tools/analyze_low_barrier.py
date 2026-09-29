"""Compare AlF barrier hypotheses using native residuals and radial wavefunctions.

Run after the recorded native calculations. This produces diagnostic artifacts;
it never selects a trial as the production model or infers a resonance lifetime.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='2'
deps=Path(__file__).resolve().parent.parent/'.pecfit_stage/plot_deps'
if deps.exists():sys.path.insert(0,str(deps))
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from coupled import fields, save_assessment
from barrier import potential, extrema, audit
from barrier_spectrum import kinetic, spectrum, inner_barrier, LIMIT, ZPE
from bound_states import read_states
from make_report import raw_rows, describe, write_csv


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('runs',type=Path);ap.add_argument('prepared',type=Path);ap.add_argument('output',type=Path)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    definitions=[('original','Selected reference','native_original','#23395b'),
                 ('above_limit_1000','Lower barrier above dissociation','refit_1000/round_04_fraction_1','#268b83'),
                 ('minimum_078','Approx. 0.78 eV above minimum','native_central_078','#d8862b'),
                 ('minimum_5800','5800 cm⁻¹ above minimum','native_exact_5800','#bd3f4c')]
    r=np.linspace(.7,8,801);T=kinetic(r);rg=np.unique(np.r_[np.linspace(.7,8,10001),np.geomspace(8,1000,1001)])
    summary={};curves={};residuals={};spectra={};reference=None
    for tag,label,relative,color in definitions:
        run=a.runs/relative;text=(run/'input.inp').read_text();fs=fields(text);field=fs['POTEN',2,2]
        if reference is None:reference=fs
        fixed={str(k):f.values==reference[k].values for k,f in fs.items() if k!=('POTEN',2,2) and k[0]!='L+'}
        if not all(fixed.values()):raise ValueError(f'Unexpected non-A curve changes: {fixed}')
        if field.values[15]!=LIMIT or field.values[24:26]!=reference['POTEN',2,2].values[24:26]:
            raise ValueError('Changed asymptote or C5/C6')
        validation=save_assessment(run,98)
        raw=raw_rows(run,a.prepared.read_text());residuals[tag]=raw
        for row in raw:row['model']=tag
        write_csv(a.output/f'{tag}_residuals.csv',raw)
        low=minimize_scalar(lambda rr:float(potential(np.array(rr),field)),bounds=(1.55,1.75),method='bounded',options={'xatol':1e-12})
        peak=inner_barrier(field,1);curve=potential(rg,field);ex=extrema(rg,curve);curves[tag]=field
        all_rows=[];counts=[];observed_caps=[]
        for j in (1,61,83,92):
            rows,b=spectrum(field,j,r,T);all_rows.extend(rows)
            counts.append(dict(J=j,below_dissociation=sum(x['below_limit'] for x in rows),
                inner_below_dissociation=sum(x['below_limit'] and x['inner_probability']>.5 for x in rows),
                inner_below_barrier=sum(x['below_barrier'] and x['inner_probability']>.5 for x in rows),
                barrier_r=b[0] if b else None,barrier_absolute_cm=b[1] if b else None))
            for v,jmax in enumerate((83,92,92,92,83,83,61)):
                if j!=jmax:continue
                inner=[x for x in rows if x['inner_v']==v]
                if len(inner)!=1:raise ValueError(f'Missing inner assignment {tag}, {j}, {v}')
                x=inner[0]
                observed_caps.append(dict(v=v,J=j,inner_probability=x['inner_probability'],energy_cm=x['energy'],
                    barrier_margin_cm=b[1]-ZPE-x['energy'],dissociation_margin_cm=LIMIT-ZPE-x['energy']))
        spectra[tag]=all_rows;write_csv(a.output/f'{tag}_spectrum.csv',all_rows)
        write_csv(a.output/f'{tag}_observed_J_caps.csv',sorted(observed_caps,key=lambda x:x['v']))
        summary[tag]=dict(label=label,source=str(run.resolve()),sha256=hashlib.sha256((run/'input.inp').read_bytes()).hexdigest(),
            minimum_absolute_cm=float(low.fun),minimum_r=float(low.x),barrier_absolute_cm=peak[1],barrier_r=peak[0],
            barrier_above_minimum_cm=peak[1]-low.fun,barrier_above_limit_cm=peak[1]-LIMIT,
            extrema=ex,finite_to_1000=bool(np.all(np.isfinite(curve))),default_shape_audit=audit(text),
            unchanged_fields=fixed,native_validation=validation,residual_statistics=describe(raw),counts=counts,
            observed_caps=sorted(observed_caps,key=lambda x:x['v']))
        config={'above_limit_1000':'shape_1000.json','minimum_078':'shape_minimum_078.json',
                'minimum_5800':'shape_minimum_5800.json'}.get(tag)
        if config:summary[tag]['relaxed_trial_shape_audit']=audit(text,json.loads((a.runs/config).read_text()))
        if tag!='original':(a.output/f'{tag}_diagnostic.inp').write_text(text)
    # Box/grid sensitivity: below-barrier inner levels must not be box artifacts.
    convergence={}
    f=curves['minimum_5800'];base=[x for x in spectra['minimum_5800'] if x['J']==1 and x['inner_v']!='' and x['below_barrier']]
    for name,rmax,n in [('grid1001',8,1001),('box12',12,1401)]:
        rr=np.linspace(.7,rmax,n);rows,_=spectrum(f,1,rr,kinetic(rr));inner=[x for x in rows if x['inner_v']!='' and x['below_barrier']]
        convergence[name]=dict(inner_below_barrier=len(inner),below_limit_finite_box=sum(x['below_limit'] for x in rows),
            max_inner_energy_change_cm=max(abs(x['energy']-y['energy']) for x,y in zip(base,inner)) if len(base)==len(inner) else None)
    native=read_states(next((a.runs/'unbound_exact_J1').glob('*.states')))
    nf=sorted([x for x in native if x['state']=='A1Pi' and x['J']==1 and x['ef']=='f' and x['energy']<LIMIT-ZPE],key=lambda x:x['energy'])
    sf=[x for x in spectra['minimum_5800'] if x['J']==1 and x['below_limit']]
    convergence['native_J1']=dict(below_limit=len(nf),flag_b=sum(x['flag']=='b' for x in nf),
        max_inner_DVR_difference_cm=max(abs(x['energy']-y['energy']) for x,y in zip(nf[:8],sf[:8])),
        max_full_DVR_difference_cm=max(abs(x['energy']-y['energy']) for x,y in zip(nf,sf)) if len(nf)==len(sf) else None)
    cv=a.runs/'native_exact_convergence'
    if (cv/'duo.out').exists():
        save_assessment(cv,98)
        cr=raw_rows(cv,a.prepared.read_text());br=residuals['minimum_5800']
        convergence['native_observed_levels']=dict(npoints=1001,vmax=180,
            max_energy_change_cm=max(abs(x['calculated_cm']-y['calculated_cm']) for x,y in zip(cr,br)))
    write_csv(a.output/'minimum_5800_native_J1_f.csv',nf)
    summary['validation']=dict(radial_rmax=8,npoints=801,inner_probability_threshold=.5,
        energy_zero='X v=0 J=0 for level tables; X potential minimum for PEC heights',
        zpe=ZPE,dissociation_absolute=LIMIT,convergence=convergence,
        count_note='Finite-box totals are lower bounds, not converged global bound-state counts. Above-barrier Pinner is not a vibrational quantum number.')
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(12,5.3),layout='constrained')
    rr=np.linspace(1.35,6,3500)
    for tag,label,_,color in definitions:
        for ax in axs:ax.plot(rr,potential(rr,curves[tag]),color=color,label=label,lw=1.7)
    for ax in axs:
        ax.axhline(LIMIT,color='.35',ls='--',lw=1,label='Common Al + F limit')
        ax.set(xlabel='Internuclear distance / Å',ylabel='Energy above X PEC minimum / cm⁻¹')
    axs[0].set(xlim=(1.4,5.5),ylim=(43000,64000),title='Lower local barriers create deep outer wells')
    axs[0].legend(fontsize=8,loc='upper right')
    axs[1].set(xlim=(1.45,2.6),ylim=(43700,52000),title='5800 cm⁻¹ trial: eight inner levels at J = 1')
    for row in base:
        axs[1].hlines(row['energy']+ZPE,1.52,1.88,color='#bd3f4c',lw=.8)
        axs[1].text(1.485,row['energy']+ZPE,str(row['inner_v']),ha='right',va='center',fontsize=8)
    axs[1].text(.97,.04,'v labels refer only to inner-well levels.\nAdditional outer/mixed bound levels remain.',
                transform=axs[1].transAxes,ha='right',va='bottom',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    fig.savefig(a.output/'barrier_comparison.png',dpi=180);fig.savefig(a.output/'barrier_comparison.pdf');plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(11,7.8),layout='constrained')
    for ax,(tag,label,_,color) in zip(axs.flat,definitions):
        rows=[x for x in residuals[tag] if x['state']==2 and x['input_weight']>0]
        sc=ax.scatter([x['J'] for x in rows],[x['residual_cm'] for x in rows],c=[x['v'] for x in rows],cmap='viridis',vmin=0,vmax=6,s=9,alpha=.65)
        st=summary[tag]['residual_statistics']['by_state']['2']
        ax.axhline(0,lw=.7,color='.4');ax.set(title=label,xlabel='J',ylabel='Observed − calculated / cm⁻¹')
        ax.text(.02,.96,f"RMS {st['rms_cm']:.3f}; median |residual| {st['median_abs_cm']:.3f}",transform=ax.transAxes,va='top',fontsize=9)
        fig.colorbar(sc,ax=ax,label='v')
    fig.savefig(a.output/'residual_comparison.png',dpi=180);fig.savefig(a.output/'residual_comparison.pdf');plt.close(fig)
    print(json.dumps({tag:{'height':summary[tag]['barrier_above_minimum_cm'],'rms_A':summary[tag]['residual_statistics']['by_state']['2']['rms_cm'],'counts_J1':summary[tag]['counts'][0]} for tag,*_ in definitions},indent=2))
    print(json.dumps(convergence,indent=2))
    table=[]
    for tag,label,*_ in definitions:
        s=summary[tag];st=s['residual_statistics']['by_state']['2'];c=s['counts'][0]
        table.append(f"| {label} | {s['barrier_above_minimum_cm']:.1f} | {s['barrier_above_limit_cm']:+.1f} | {c['inner_below_dissociation']} | {st['rms_cm']:.4f} | {st['median_abs_cm']:.4f} |")
    cv=convergence.get('native_observed_levels',{}).get('max_energy_change_cm')
    report=f'''# AlF: testing a lower A-state barrier

The requested hypothesis was tested, including direct native Duo refits and constrained A-state refinements. A smooth trial with eight inner-well levels (v = 0–7 at J = 1) was obtained, but its experimental agreement is much worse and it creates a deep outer well containing additional bound states. **It is a diagnostic counterexample, not an accepted replacement model.** The selected X/A model is unchanged.

## Fit and state-count comparison

All heights below are measured from the actual adiabatic A minimum or the common Al + F dissociation limit, as explicitly labelled. The nominal 1000 cm⁻¹-above-limit trial relaxed to 894.2 cm⁻¹ in its allowed fitting window. The approximate 0.78-eV trial reached 6310.8 cm⁻¹ (0.78244 eV) above the actual minimum. The final eight-level trial reached 5800.17 cm⁻¹ (about 0.7191 eV).

| Model | Barrier above A minimum / cm⁻¹ | Barrier above dissociation / cm⁻¹ | Inner levels below dissociation, J=1 | A RMS / cm⁻¹ | A median absolute residual / cm⁻¹ |
|---|---:|---:|---:|---:|---:|
{chr(10).join(table)}

The table uses all 1097 positive-weight original A records, retaining the six repeated assignments and excluding only the original twelve zero-weight records from the statistics. All 1109 original A records and 731 X records are supplied in each residual file. No additional observations were discarded. Native fits use 1103 unique A levels, including the zero-weight levels; duplicate means are expanded back to the source records for this comparison. Every native evaluation matched all supplied assignments. X PEC, X BOB, the absolute limit and A C5/C6 were kept fixed. X RMS remains approximately 0.0027 cm⁻¹; its energies can shift slightly through X/A coupling.

![Barrier alternatives](barrier_comparison.png)

![All-data residual comparison](residual_comparison.png)

## What “vmax = 7” means in this trial

The final barrier is at r = 2.06085 Å and 49762.48 cm⁻¹ above the X potential minimum. The common dissociation limit is 55564.02471 cm⁻¹. Consequently the curve must rise again after the local maximum to reach its asymptote. With the retained attractive C5/C6 tail, the lower-barrier construction develops an outer minimum at approximately 2.241 Å and 49404.10 cm⁻¹ — about 6160 cm⁻¹ below dissociation. This is a substantial change in the physical PEC, not a shallow van der Waals well.

At J = 1 the final model has eight levels with more than half their probability inside the local barrier and energy below that barrier. Their inner probabilities range from 0.9685 to essentially 1. There are nevertheless **44 A(f) eigenlevels below dissociation in the 8 Å box**, and 45 in a 12 Å box. The latter is a lower bound on the full bound spectrum, not a claim that the diffuse outer spectrum is converged. Above the local barrier, “inner v” based on a 50% probability cutoff ceases to be a reliable spectroscopic assignment.

Native Duo with the requested FITTING-off INTENSITY/UNBOUND/STATES_ONLY block confirms the 44 below-limit f-parity levels at J = 1. Of those, **37 carry the b flag**, including outer/mixed levels; only eight are the inner-well vibrational progression. Thus THRESH_BOUND=0.01 and THRESH_bound_rmax=4 cannot enforce a maximum inner v by themselves. The printed v label is also not a general inner-well quantum number. Counts here exclude M and nuclear-spin degeneracies and refer to one f parity unless stated otherwise.

## Observed rotation and numerical checks

The source reaches A v = 6 and J = 92 overall; v = 6 reaches J = 61. At each observed per-v maximum J, the final trial still has an inner-localized level below the rotational barrier. For v=6, J=61, its inner probability is 0.9389 and it is 340.3 cm⁻¹ below that effective barrier. For v=5, J=83, the values are 0.9946 and 512.5 cm⁻¹. These are localization checks of model levels, not successful fits: the corresponding energy residuals remain large. See the per-model `*_observed_J_caps.csv` tables.

The A(f) radial calculation uses the full uncontracted sinc DVR with 801 points, the native AlF reduced mass, and Duo's J(J+1)−2 convention for this supplied singlet-Pi Hamiltonian. Full X/A native Duo handles both parities for every observed level through J=98. Increasing the final trial to 1001 radial points and vmax=180 changes the observed energies by at most {cv:.3g} cm⁻¹. Therefore the roughly 18.57 cm⁻¹ A RMS is not a basis-convergence error.

The eight low-J inner levels survive both the 1001-point check and extension to a 12 Å, 1401-point box. Native six-decimal `.states` energies for those levels agree with the independent radial calculation to {convergence['native_J1']['max_inner_DVR_difference_cm']:.3g} cm⁻¹. Diffuse outer levels are less converged: across the full 44-level set, the 801/1001 discretizations differ by up to {convergence['native_J1']['max_full_DVR_difference_cm']:.4g} cm⁻¹, and the larger box adds another level. No radiative lifetime, tunnelling width or predissociation lifetime has been inferred.

An earlier 5800-cm⁻¹ trial also received the complete native J=0–98 UNBOUND calculation (53,090 states). Its A RMS was 19.34 cm⁻¹; all 1103 supplied unique A assignments were localized and below the common limit. This is retained as separate evidence and is not confused with the final minimum-referenced trial, for which the native UNBOUND check is J≤1 and the full experimental energy validation is J≤98.

## The 0.763–0.78 eV reference

[Andreazza & de Almeida (2014)](https://academic.oup.com/mnras/article/437/3/2932/1034811) report a 0.763 eV hump at 4.8 bohr and discuss it as a barrier suppressing association of incoming Al + F atoms at low collision energy. Their collision formula uses the potential relative to separated atoms. This supports interpreting the reported hump above dissociation, rather than above Te. The sentence giving the hump does not explicitly restate its zero; the interpretation is inferred from the collision-energy context, not a newly recovered quotation from the original 1988 paper.

[Qin, Bai & Liu (2022), Table 1](https://academic.oup.com/mnras/article/510/2/3011/6460490) independently assign X and A to the same ground-atom asymptote. A 0.78-eV barrier above the present A minimum would be near 50241 cm⁻¹ on the X-minimum scale, already about 5323 cm⁻¹ below that common dissociation limit. The retained X/A separation to dissociation cannot be removed merely by lowering a local barrier. The earlier literature note has been clarified to distinguish the contextual energy-zero inference from the reported number.

The highest observed v is an observational limit; it does not establish that higher vibrational levels are absent. A fit constrained to seven or eight *total* bound levels would require additional, independently justified changes to the potential and/or dissociation physics. The numerical trials here do not prove that every possible alternative functional form fails; they show that the tested low-barrier versions of this two-state, fixed-asymptote model are unsuitable replacements.

## Methods, files and reproduction

- The lower-above-dissociation trial used four native robust-fitting rounds, three proposal iterations each, ROBUST=0.00001, with fresh all-J validation and an explicit relaxed shape profile. Its robust median improves while its ordinary RMS increases slightly; it still has 17 inner bound levels.
- The below-dissociation trials began by matching the inner PEC while imposing a lower peak, then fitting experimental energies. Ordinary native proposals lost the required local maximum and were rejected. A constrained A(f) Hamiltonian fit then varied A Te/Re/B0–B6 and, in additional trials, repulsive amplitude and crossing parameters. Equal PL=PR=6 and fixed C5/C6 were retained. Native e/f corrections were frozen only inside that surrogate; every reported residual uses a fresh complete native Hamiltonian. The final eight-level trial constrained the actual adiabatic minimum, not the auxiliary EMO Te or Re value. The approximate 0.78-eV run hit its iteration limit and is labelled exploratory, not a proven optimum.
- Shape profiles allowing outer wells several thousand cm⁻¹ deep were used solely to explore the hypothesis. Passing those relaxed profiles does not certify physical adequacy. All published trial PECs remain finite on 0.7–1000 Å; the original production shape constraints remain the default and reject these alternatives.
- `minimum_5800_diagnostic.inp` is the final eight-inner-level trial; `minimum_078_diagnostic.inp` and `above_limit_1000_diagnostic.inp` are comparison hypotheses. These are complete native Duo inputs, clearly separated from `refined_model/27AlF_X_A_coupled_pec.inp`.
- `summary.json`, spectra, original-record residuals, observed-J-cap tables, PNG/PDF figures and native J=1 f-level tables accompany this report. Full inputs, stdout, .states, manifests, trial histories and optimizer snapshots are retained under `AlF/pecfit_runs/low_barrier/`.
- Reusable tools added: `low_barrier_seed.py`, `fit_low_barrier.py`, `analyze_low_barrier.py`; `refine.py --shape-config` accepts explicit trial bounds without changing production defaults. `barrier_spectrum.py` now scans at least to dissociation even if the local barrier is lower. The native localization summary no longer calls all b-labelled states inner-well states. Ten regression checks pass.

For a fresh native validation, use a new output directory:

```powershell
python AlF/refine_tools/run_native.py AlF/barrier_trials/minimum_5800_diagnostic.inp AlF/pecfit_runs/low_barrier_recheck --duo C:/sergei/programs/Duo/build/ifx/duo.exe --jmax 98 --precise
python AlF/refine_tools/coupled.py AlF/pecfit_runs/low_barrier_recheck --jmax 98
python AlF/refine_tools/bound_states.py AlF/barrier_trials/minimum_5800_diagnostic.inp AlF/pecfit_runs/low_barrier_unbound_recheck --duo C:/sergei/programs/Duo/build/ifx/duo.exe --jmax 1 --npoints 1001 --vmax 180
python AlF/refine_tools/barrier_spectrum.py AlF/barrier_trials/minimum_5800_diagnostic.inp AlF/pecfit_runs/low_barrier_radial_recheck --js 1 61 83 92 --npoints 801
```

The recommended next constraint is independent near-dissociation or lifetime evidence, rather than treating v=6 as an enforced dissociation cutoff. The current selected model remains available unchanged while those physical constraints are resolved.
'''
    (a.output/'LOW_BARRIER_REPORT.md').write_text(report,encoding='utf-8')


if __name__=='__main__':main()
