# Diagnostic fitting step with the Qin reference

This is a bounded diagnostic of the newly added constraint, **not a converged
refinement and not evidence that an acceptable joint fit is impossible**.

The preferred input uses the 66 supplied A-state points, fixed initial
PS1997 weights (alpha=0.001, cutoff=18000 cm-1) and fit_factor=100. All synthetic
extrapolation points have zero weight. The prepared inputs retain robust=1e-5;
for this diagnostic only, robust weighting was disabled to measure the
ordinary least-squares response to the explicit constraint. The native update
scale was 0.1. The first complete proposed parameter block was extracted and
evaluated in a separate zero-iteration run through J=98.

## A-only trial

Only the nine already marked A-PEC parameters were free: the auxiliary EMO
V0, RE and B0--B6. X, the electronic coupling, BOB, the repulsive amplitude,
coupling-function construction parameters and the common dissociation limit
were held at their previous values. No new observation was discarded.

| Metric / cm-1 | Baseline | First A-only proposal |
|---|---:|---:|
| A experimental weighted RMS | 0.582176 | 25.939067 |
| X experimental weighted RMS | 0.002672 | 0.002676 |
| Combined experimental weighted RMS | 0.450996 | 20.094147 |
| Reference PEC weighted RMS | 2604.635136 | 2416.728351 |

The experimental metrics use the same original input weights on the same
1822 positive-weight prepared observations (731 X and 1091 A); duplicates
were already combined in that prepared input. These numbers must not be
confused with an expanded residual statistic on every original source record.
Reference RMS uses the fixed PS1997 point weights, normalized within the
reference set. The two kinds of RMS are not interchangeable tolerances.

The A-only candidate passes the existing shape guard but is **not selected**:
its experimental fit is much worse. Its barrier is about 60481.7 cm-1 at
2.5717 A, still far above the Te-aligned source curve. Holding the old outer
construction parameters and dissociation limit fixed substantially restricts
how the fit can follow the supplied reference. The existing shape guard also
retains the older above-asymptote barrier range; passing it does not demonstrate
agreement with Qin's supplied barrier or establish physical accuracy.

An earlier unrestricted diagnostic, allowing all inherited fit flags, also
reduced reference RMS but distorted X and the corrections and failed the
shape guard. Its combined experimental RMS rose to about 1338.94 cm-1. That
proposal was rejected; its raw output is retained solely for audit.

## Consequence

The reference is active and parsed correctly. A useful continued refinement
needs an explicitly balanced joint objective, controlled parameter bounds,
and a decision about the source grid's well-depth discrepancy and common
dissociation limit. Blindly adding fit_factor=100 to the previous unconstrained
parameter set is not a validated route to higher spectroscopic accuracy.
No inference about maximum bound v/J or lifetimes follows from these trials.

## Evidence and reproduction

`trial_result.json` contains the A-only shape audit and fresh native residual
statistics. Complete raw runs are under `AlF/pecfit_runs/qin_2022/`, including
the native `fit.pot`, input/executable hashes, full-precision parameter
proposals and independent `rovibronic_energies.dat` evaluations.

From the repository root, use a new directory:

```powershell
python AlF/refine_tools/trial_abinitio_reference.py AlF/abinitio/Qin_2022/inputs/AlF_with_reference_fixed_weights.inp AlF/pecfit_runs/qin_trial_new --duo C:/sergei/programs/Duo/build/ifx/duo.exe --scale 0.1
```

The native `itmax=1` run emits two complete proposal blocks. The diagnostic
deliberately assesses block zero, the first update, rather than mistaking the
residual table for evidence about an unevaluated final parameter proposal.
