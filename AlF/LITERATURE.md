# AlF physical assumptions and literature

The model has exactly two electronic states: X1Sigma+ (Lambda = 0, multiplicity 1) and A1Pi (Lambda = 1, multiplicity 1). The CH example supplies the **functional form**, not AlF numerical parameters or symmetry.

## Fixed long-range information

[Qin, Bai & Liu, MNRAS 510, 3011 (2022), Table 1 and equations 7–9](https://doi.org/10.1093/mnras/stab3598) assign X and A to the same ground-atom Al(2P) + F(2P) limit. They estimate C6 = 75.23 with the London formula and fit C5 = 8.07 for A. These are approximate theoretical estimates, not precision spectroscopic determinations. The table does not explicitly label coefficient units: the atomic-unit interpretation of C6 follows their London formula and stated atomic-unit polarizabilities; C5 is interpreted in the consistent atomic-unit convention.

Using Eh = 219474.63136320 cm-1 and a0 = 0.529177210903 Angstrom gives:

| Coefficient in the AlF input | Value | Units |
|---|---:|---|
| B5 = -C5 | -73495.9702828864 | cm-1 Angstrom^5 |
| B6 = -C6 | -362561.9198027829 | cm-1 Angstrom^6 |

The supplied absolute dissociation energy, 55564.02471 cm-1 relative to the X potential minimum, is retained. It agrees with the fixed AlF value in [Yousefi & Bernath, ApJS 237, 8 (2018)](https://doi.org/10.3847/1538-4365/aacc6a). It is not refitted from the low vibrational levels.

## Barrier and coupling

[Andreazza & de Almeida, MNRAS 437, 2932 (2014)](https://academic.oup.com/mnras/article/437/3/2932/1034811) report an A-state hump of 0.763 eV at 4.8 bohr, using earlier electronic-structure curves. The energy converts to approximately 6154.01 cm-1 and the radius to 2.54005 Angstrom. Their discussion of suppressed association at low collision energy, and their collision formula referenced to separated atoms, support interpreting the hump above dissociation. This energy-zero assignment is a contextual inference: the sentence reporting the hump does not explicitly restate its zero, and the original 1988 calculation has not been independently recovered here. These are initialization targets, not additional measured levels.

The [lower-barrier tests](barrier_trials/LOW_BARRIER_REPORT.md) also examine the user's alternative interpretation above the A minimum. With the retained common Al + F asymptote, those trials introduce a deep outer well. One supports only v=0-7 in its inner well at J=1, but many additional outer/mixed bound levels remain and the experimental fit deteriorates sharply. The original selected model has not been replaced on this basis.

The native Duo implementation used here is:

```
Vrep(r) = D + A exp(-delta/r)/r^gamma + sum_k Bk/r^k
VA(r) = (VEMO(r) + Vrep(r) - sqrt((VEMO(r)-Vrep(r))^2 + 4 W(r)^2))/2
```

W is a second EMO with zero asymptote. `COMPON 1` selects the lower branch. `Nparameters` describes sub-functions within **one** electronic PEC; it does not add electronic states to the Hamiltonian.

Gamma is fixed at 8, so the positive repulsive term decays faster than the retained r^-5 and r^-6 attractions. Copying CH's gamma = 5.5 would make that term dominate C6 at sufficiently large distance. The diabatic EMO asymptote (80000 cm-1), repulsive amplitude, delta, and W parameters are construction choices constrained by the barrier target and smoothness. They are not independently measured AlF constants. The shallow outer van der Waals well follows from a barrier followed by an attractive tail.

## Scope of the empirical fit

The observed levels sample X v = 0–9 and A v = 0–6; they cannot uniquely determine the barrier or long-range tail. PL and PR are equal to 6 in all EMO sub-functions. X remains an EMO: its limiting energy is fixed, but its exponential tail is not a dispersion expansion.

The interaction of A1Pi with b3Sigma+ is documented by [Doppelbauer et al., Molecular Physics (2020)](https://doi.org/10.1080/00268976.2020.1810351). A two-state X/A model cannot explicitly represent those triplet perturbations. Large or localized residuals are flagged for data/perturbation review; their cause is not assigned automatically, and they remain in reported all-data statistics.

## Supplied Qin 2022 numerical reference

The complete [Qin reference package](abinitio/Qin_2022/README.md) preserves the
user-supplied source file and all 66 A1Pi radius/energy pairs. Each state has
its own radial column. The A grid is already in Angstrom; energies are treated
as cm-1. Its sampled minimum is aligned to the selected model's actual A
minimum with a constant shift of +53721.79358272087 cm-1.

The supplied grid has a sampled barrier 10587.52 cm-1 above its minimum and
842.14 cm-1 above its 9-Angstrom endpoint. Its endpoint-based well depth,
9745.38 cm-1, differs substantially from De(A)=13682.94 cm-1 in Qin's Table 2.
Its sampled A/X minimum difference, 47917.98 cm-1, also differs from the
table's Te=44542.07 cm-1. This discrepancy is unresolved; the source grid has
not been silently rescaled to match the table or the retained empirical limit.

After Te alignment the endpoint is 1868.65 cm-1 below the model's fixed common
dissociation limit. Native syntax and weights have been checked through J=98.
The [limited fitting trials](abinitio/Qin_2022/TRIAL.md) show the reference is
active but do not improve the accepted fit. The selected model remains unchanged.
The prepared explicit-weight variant avoids this Duo build assigning positive
PS1997 weights to synthetic short-range extrapolation points.

## Uniformly scaled Qin A-state test

The [scaled-Qin report](abinitio/Qin_2022_scaled/README.md) tests scaling the
supplied numerical A well to the adopted common dissociation limit. The
file-derived depth gives a factor 1.18852526; the finite-endpoint convention
gives 1.19174738. Both native GRID tests retain 17 inner bound levels at J=1.
A weak full-curve constraint lowers the analytic barrier and approximately
preserves experimental accuracy, but does not reduce vmax to 6 or 7. The
file/table De discrepancy and the uncertainty of the outer diffuse spectrum
remain explicit. The selected model is unchanged.
