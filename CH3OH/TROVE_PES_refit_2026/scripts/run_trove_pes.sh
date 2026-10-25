#!/bin/bash
# run_trove_pes.sh <params.trove> <rundir> [jmax]
# Exact potential.chk + consistent EQUILIBRIUM/SPECPARAM, then TROVE J=0..jmax.
set -e
TROVEFILE=$(readlink -f $1); RUN=$2; JMAX=${3:-3}
A=/home/user/ExoMol/CH3OH/pes_refit/analysis
EXE=/home/user/ExoMol/CH3OH/pes_refit/trove/src/j-trove.x
mkdir -p $RUN && cd $RUN
cp /home/user/ExoMol/CH3OH/nb/{kinetic.chk,combination.chk} .
python3 $A/generate_ch3oh_potential_chk.py $TROVEFILE --out potential.chk --thresh 1e-8 | tee potential_chk.log
python3 $A/make_trove_inputs.py --coords potential.chk.trove_coords --poten $TROVEFILE --jmax $JMAX
for J in $(seq 0 $JMAX); do
  /usr/bin/time -f "J=$J: %e s" env OMP_NUM_THREADS=32 $EXE < J$J.inp > J$J.out 2>&1 || true
  tail -1 J$J.out
  python3 $A/trove_energies.py J$J.out | tail -4
done
