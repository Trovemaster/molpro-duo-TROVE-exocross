# AlF refinement tools

This driver builds a two-state X1Sigma+/A1Pi model. The A potential uses `coupled-pec`, with sub-types `EMO repulsive_exp EMO` and the lower adiabatic component. The AlF-specific long-range estimates, barrier target, and their limitations are in `../LITERATURE.md`.

Install the scientific dependencies in your Python environment:

```powershell
python -m pip install -r AlF/refine_tools/requirements.txt
```

Run the complete workflow from the repository root, using a new output directory:

```powershell
python AlF/refine_tools/pipeline.py AlF/reference_inputs/27AlF_X1E+_A1Pi_EMO_02.inp --duo C:/sergei/programs/Duo/build/ifx/duo.exe --output AlF/pecfit_runs/new_fit
```

The workflow audits the source, constructs the barrier seed, fits to experiment, checks radial-grid/basis/box convergence, evaluates every retained original observation, and writes `final.inp`, plots, and JSON diagnostics. A full run takes several minutes or longer depending on hardware and the number of refinement rounds. `--ordinary-rounds` and `--robust-rounds` control the outer loops. Publication-quality residuals and final parameters must be reviewed together with the flagged experimental levels.

For a fresh evaluation of the delivered model:

```powershell
python AlF/refine_tools/run_native.py AlF/refined_model/27AlF_X_A_coupled_pec.inp AlF/pecfit_runs/check --duo C:/sergei/programs/Duo/build/ifx/duo.exe --jmax 98 --precise
python AlF/refine_tools/coupled.py AlF/pecfit_runs/check --jmax 98
python AlF/refine_tools/check_curves.py AlF/pecfit_runs/check
```

To refine further from an existing model:

```powershell
python AlF/refine_tools/refine.py AlF/refined_model/27AlF_X_A_coupled_pec.inp --duo C:/sergei/programs/Duo/build/ifx/duo.exe --output AlF/pecfit_runs/continued_fit --jmax 98 --rounds 4 --proposal-steps 3 --scale 1 --robust 0.00001 --criterion cauchy
```

The delivered grid is larger than the initial fitting grid, so this command may run more slowly. Every proposed parameter set is checked for sensible PEC/BOB/coupling shapes and then recalculated with native Duo before acceptance. `weighted_rms` is the alternative acceptance criterion for ordinary least squares. Cauchy acceptance may increase ordinary RMS while improving the robust objective; both are retained in the history.

## Data handling

The original file is retained unchanged. `input_audit.json` records character normalization, the missing Lx morphing reference, the singlet-Pi Omega correction, the two source records already disabled by negative J, the extension of JLIST to 98, and six conflicting duplicate assignments. Duplicate levels are represented by their weighted mean and summed weight during fitting. Final residuals expand these means back to all original records. The original twelve zero-weight levels remain zero-weight. No additional observation is rejected because its residual is large.

The duplicates may contain parity transcription errors. The preparation deliberately does not guess different parities. A robust loss on the duplicate means is not mathematically identical to a robust loss on separate duplicate measurements; the ordinary weighted squared loss differs only by a fixed within-group constant.

## Safeguards and reproducibility

- Compound parameter labels can repeat. Parameter transfer uses object, state indices, and position; it never identifies a coefficient by `B6` alone.
- Native parameter tables contain *proposed* values. Their preceding residual table is not evidence for those values. Each accepted proposal receives a fresh zero-iteration evaluation.
- Restarted robust fits carry Watson IRLS weights explicitly. Native Duo updates robust weights after its first step, so restarting one-step fits with only `ROBUST` would otherwise repeat ordinary least squares.
- Every run records input/executable hashes and preserves stdout. The refinement also snapshots its Python source.
- The selected PECs are checked from 0.7 to 1000 Angstrom. X must have one well; A must have one inner well and one barrier, with at most one shallow outer well. Equal PL/PR = 6 or 8 is enforced. Small rotational and electronic angular-momentum corrections are bounded.
- `convergence.py` compares 601/801 radial points, vmax 60/90, and a final 1001-point, 8-Angstrom, vmax-100 calculation on the same levels. Its pass tolerance is 0.0001 cm-1 per positive-weight level.
- Scientific tests: `python -m unittest discover -s AlF/refine_tools -p test_refinement.py -v`.

The prepared inputs use integer J, numeric state labels, and the supplied AlF level conventions. The driver, literature constants, and shape bounds are specific to AlF; revise them explicitly before applying this barrier model to a different molecule. NumPy, SciPy, Matplotlib, and a Duo executable that supports `coupled-pec`/`repulsive_exp` are required. The input has no intensity calculation, so its inherited dipole functions are not validated by this energy fit.

## Bound-state and barrier diagnostics

`bound_states.py` removes FITTING and runs the requested INTENSITY / UNBOUND /
STATES_ONLY calculation. It records native `.states`, mean radii and tail
probabilities, subtracts the actual Duo ZPE from the coupled PEC asymptote, and
checks the supplied observations. It preserves PEC and coupling parameters.
The diagnostic distinguishes energy below dissociation from Duo's b/u
localization flag and warns when bound levels still reach the requested J cap.

```powershell
python AlF/refine_tools/bound_states.py AlF/refined_model/27AlF_X_A_coupled_pec.inp AlF/pecfit_runs/new_bound_scan --duo C:/sergei/programs/Duo/build/ifx/duo.exe --jmax 160 --vmax 180
```

Append `--bound-diagnostics` to `pipeline.py` to run this stage automatically
after refinement; `--bound-jmax` defaults to 160 and `--bound-vmax` to 180.
The native driver defaults to the user's J <= 40, 8 A box, 1001 radial points,
THRESH_BOUND=0.01 and THRESH_bound_rmax=4 when run directly. Larger boxes need
both adequate radial resolution (`--npoints`) and enough contracted functions
(`--vmax`); the above-threshold spectrum is not converged merely because the
fitted low-lying levels are converged. Use a new output directory for each run.

`barrier_spectrum.py`, `check_outer_well.py` and `make_bound_report.py` reproduce
the additional AlF-specific radial cross-checks and plots. Their single-channel
A(f) calculation uses the present input's mass, diagonal-L2 convention and
energy zero; it is not a substitute for general multistate Duo or a resonance
width solver. Consult [the report](../bound_diagnostics/BOUND_STATE_REPORT.md)
before interpreting high printed v labels or outer-well levels. These routines
do not calculate radiative or predissociation lifetimes.

Regression checks: `python -m unittest discover -s AlF/refine_tools -p test_bound_states.py -v`.

## Lower-barrier sensitivity trials

See [the completed comparison](../barrier_trials/LOW_BARRIER_REPORT.md). The
trial inputs in `AlF/barrier_trials/` are diagnostic hypotheses, not accepted
replacements for the selected model. Raw runs are in `AlF/pecfit_runs/low_barrier/`.

`refine.py --shape-config PATH.json` overrides only the named shape bounds for
an explicitly documented trial. Unknown keys are rejected; production defaults
remain unchanged. The published profiles deliberately allow deeper outer wells
to test the hypothesis, so passing them does not establish physical adequacy.

`low_barrier_seed.py INPUT OUTPUT --height 5800 --reference minimum` creates
an initial lower-barrier curve. It does not fit the experimental energies.
`fit_low_barrier.py BASELINE_DIR OUTPUT --height 5800 --minimum-reference
--flexible-crossing --loss-scale 1` then fits the experimental A levels with
the local barrier inside the objective. BASELINE_DIR must contain a complete
native `input.inp` and its `residuals.csv` from `coupled.py`. This experimental
surrogate is restricted to the current AlF 15/12/13 parameter layout and fixed
X PEC. It freezes native e/f corrections internally: every proposed model
requires fresh full-J native validation and radial/box checks. It never
automatically promotes a result to the production model. Use `--minimum-reference`
to measure the height from the actual adiabatic minimum; the legacy default
reference is V(RE_EMO), which need not equal that minimum. `--jacobian relative`
replays the early forward-difference trials; central differences are the default.

`analyze_low_barrier.py AlF/pecfit_runs/low_barrier
AlF/pecfit_runs/low_barrier/prepared_original_records.inp OUTPUT` reproduces
the native-residual comparison and plots from the saved trials.

The independent radial spectrum always extends at least to dissociation,
even when a local barrier lies below it. `bound_states.py` reports
`J_cutoff_reached_by_localized_below_limit_levels`; the earlier field name
incorrectly implied that all localized levels were inner-well levels.
Outer-well and mixed states can satisfy Duo's b flag and mean-radius cutoff.
Use wavefunction probability inside the local barrier to identify the inner
progression, and compare larger radial boxes for diffuse bound levels.

## Ab initio reference preparation

The [Qin 2022 reference package](../abinitio/Qin_2022/README.md) contains the
complete supplied A1Pi grid, Te-aligned Duo inputs, comparison plot, source
hashes and native validation. `prepare_qin_reference.py` extracts the correct
independent A radius/energy pair and calls `add_abinitio_reference.py` for the
unit/energy-reference conversion. It requires an explicit energy unit.

`trial_abinitio_reference.py` performs one A-only diagnostic update and a
fresh full-J assessment without promoting a model. Read the package's source
discrepancy and weighting notes before fitting. The existing `refine.py`
acceptance objective is experimental-only and has not been extended here to
a joint experimental/reference objective.

Data-integrity tests: `python -m unittest discover -s AlF/refine_tools -p test_abinitio_reference.py -v`.

## Uniform scaling and grid-state diagnostics

See [Qin_2022_scaled](../abinitio/Qin_2022_scaled/README.md) for the completed
native grid, constrained-refinement and convergence experiments.
`scale_qin_reference.py` prepares both energy-zero conventions;
`bound_states.py --dissociation` supports a GRID A PEC with an explicit
absolute asymptote; `grid_bound_spectrum.py` resolves inner/outer character
from native potential grids. `fit_scaled_qin.py` implements a joint
spectroscopic/reference surrogate objective followed by mandatory native
validation. Its reference weights and strengths are documented explicitly;
it does not replace the experimental-only objective in `refine.py`.
The current selected model is preserved. Tests include energy-zero, scaling,
reference-template isolation and bound/localization distinctions.
