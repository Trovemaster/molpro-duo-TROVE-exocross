# CH3OH: refitted ab initio PES and TROVE calculations (October 2026)

This folder holds a refit of the CH3OH potential energy surface (PES) and the TROVE calculations
built on it:

* **PES.** Fitted to 49,537 CCSD(T)/CBS energies (`CH3OH-CBS`).
* **KEO.** The exact kinetic energy operator from `AB3CD.nb`.
* **Coordinates and symmetry.** TROVE `R-ALPHA132-THETA-S-TAU` coordinates and the `C3v(M)-2`
  representation.

All TROVE runs (J = 0–3, test bases up to production-type bases) finish **without
symmetry errors**.

## Layout

| Path | Content |
| --- | --- |
| `code/CH3OH_fitpes.f90` | Fitting program. It follows the structure of the CH3Br fitting program and uses the CH3OH coordinates and the C3v(M)-2 operations. The input format is the same as for `CH3OH_v10_theta12_14-MAPLESym-13-09.f90`. `build.sh` builds it. `pes_eval.f90` is a fast OpenMP PES evaluator. |
| `fit/CH3OH_PES_fitW40_r2.trove` | **Final PES** in TROVE `POTEN` format (`poten_ch3oh_sym_IV` term convention), with the matching `EQUILIBRIUM`/`SPECPARAM` blocks. |
| `fit/CH3OH_fitW40_r2.inp.gz` | Fit input: 8 fixed non-linear parameters, 4,839 terms, the ab initio data with weights, and 35 hole pseudo-points. |
| `fit/CH3OH_fitW40_r2.out.gz`, `fit/00.res.gz` | Fit log and residuals. |
| `scripts/` | Python tools: coordinates and symmetry (`ch3oh_pes.py`), term generation, fit drivers, hole checks, the exact `potential.chk` generator (`generate_ch3oh_potential_chk.py`), the TROVE input generator and the analysis tables. The scripts contain absolute paths of the original working tree (`/home/user/ExoMol/CH3OH/...`); edit them before re-running. |
| `trove/potential.chk`, `trove/potential.chk.trove_coords` | TROVE `potential.chk` ("with_modes") for the final PES. It is an exact analytic conversion: 23,616 terms with a maximum error of 3×10⁻¹⁰ cm⁻¹. The `.trove_coords` file holds the matching `EQUILIBRIUM`/`SPECPARAM`. |
| `trove/kinetic.chk`, `trove/combination.chk` | Exact KEO from `AB3CD.nb`, as used in all runs. |
| `trove/test_frozen_bends/` | J = 0–3 test runs with the bends frozen (stretches v ≤ 2). `old_pes_consistent_J0/` re-runs the old PES with its own fitted `EQUILIBRIUM`/`SPECPARAM` for comparison. |
| `trove/bends_unfrozen/` | J = 0–3 with the bends unfrozen: P6, P8, P6T8, P6B4, plus a J = 0-contraction validation. |
| `trove/production/` | Production-type bases: PB6, PB6T11 and PB8E10 (40,214 functions). J = 1–3 are run in the J = 0 contracted basis using the sliced workflow. |
| `trove/cluster_pipeline/` | Inputs for running the J = 0-contraction pipeline as consecutive batch jobs in one working directory: `make_pipeline_inputs.py`, the P6 pipeline test and the PB8T11 production basis. `BUILD_INFO.md` describes the TROVE build. |
| `trove/tools/` | Basis generators (`make_prod_inputs.py`, `make_slice_inputs.py`, `make_j0contr_inputs.py`), local runners and basis-size estimators. Also the OpenBLAS shim `blas_shim4.c`, which is only needed on machines without a system OpenBLAS. |
| `docs/REPORT_pes_refit_trove_production.md` | Project report: fit statistics, hole checks, band origins and the production-basis convergence study. |
| `figures/` | 1D cuts of the PES and the fit residuals. |

## Important conventions

1. **`EQUILIBRIUM`/`SPECPARAM` must match the fit.** TROVE builds the `morse(x-xe)` and `x-xe`
   potential coordinates of `potential.chk` from its `EQUILIBRIUM` and `SPECPARAM` blocks, so these
   must be exactly the fitted values. They are given in the `.trove` file and in
   `potential.chk.trove_coords`.

   The working input used before this work broke this rule. As a result, the PES that TROVE used
   was distorted by 103 cm⁻¹ rms below 5000 cm⁻¹, and the CH stretches came out about 150 cm⁻¹ too
   low.
2. **Torsional basis weight.** Use `resc 0.34` for the torsional (Fourier) basis. This gives one
   polyad unit per torsional triplet; `range 0,8` gives 9 functions and `range 0,11` gives 12.
3. **J = 0 contraction.**
   * `ERANGE` must come after `ENERCUT` in the `CONTRACTION` block, because `ENERCUT` also sets
     `erange(2)`.
   * J > 0 runs need `ZPE <value>` in the `DIAGONALIZER` block. Take the value from the J = 0
     output (`Zero-point-energy is ...`).
4. **TROVE version.** Use Trovemaster/TROVE `develop` at commit 930f7f6 or later. It contains the
   CH3OH symmetry fixes 70e668d and 1ecc400, and the fix 369ef18 for the J = 0 crash in
   `PThamiltonian_contract`.

## Re-running

```bash
# Fit the PES
cd code && ./build.sh
zcat ../fit/CH3OH_fitW40_r2.inp.gz > fit.inp
OMP_NUM_THREADS=48 ./CH3OH_fitpes_ob.x < fit.inp > fit.out   # writes 00.res and CH3OH_PES_fitW40_r2.trove

# Exact potential.chk for TROVE (also writes potential.chk.trove_coords)
python3 ../scripts/generate_ch3oh_potential_chk.py CH3OH_PES_fitW40_r2.trove --out potential.chk

# Production-type TROVE inputs (v <= 4 in all classes, torsion 0..11, polyads <= 8)
python3 ../trove/tools/make_prod_inputs.py ../trove/bends_unfrozen/P6/J0.inp RUN --pmax 8 --tors 11
```

### J = 0-contraction pipeline in one working directory

`trove/cluster_pipeline/make_pipeline_inputs.py RUN/J0.inp RUN PB8T11 --erange 8000` writes three
inputs:

* `PB8T11_J0s.inp`: the J = 0 energies, the vibrational matrix elements and the eigenvectors.
* `PB8T11_slices.inp`: the 12 rotational/Coriolis matrix-element slices, `matelem1..12.chk`.
  Each slice is N × N × 8 bytes, about 49 GB for the 78,234 functions of PB8T11.
* `PB8T11_J0conv.inp`: the conversion to the J = 0 eigenbasis.

After `J0s`, run the script again with `--zpe <ZPE> --jmax J` to get the `J<n>c.inp` inputs for
J > 0.

Resources for PB8T11 J0s: about 100 GB of RAM in use (2 N² × 8 bytes) and 48 cores; about 50 GB of
local scratch is needed for `contr_matelem.chk`.

The pipeline was checked on the P6 test basis. A cluster build of TROVE (see `BUILD_INFO.md`)
reproduces all 4,222 J = 0 energies and the ZPE (10893.391847 cm⁻¹) of the local build. The full
chain (J0s, slices, J0conv, J1c) runs without symmetry errors.
