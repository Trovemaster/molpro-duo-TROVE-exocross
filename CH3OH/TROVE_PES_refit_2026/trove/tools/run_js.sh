#!/bin/bash
# run_js.sh <rundir> <jmin> <jmax> [threads]: run prepared TROVE inputs J<jmin>..J<jmax>.inp
RUN=$1; J1=$2; J2=$3; NT=${4:-32}
A=/home/user/ExoMol/CH3OH/pes_refit/analysis
EXE=/home/user/ExoMol/CH3OH/pes_refit/trove/src/j-trove.x
cd $RUN
for J in $(seq $J1 $J2); do
  /usr/bin/time -f "J=$J: %e s, maxRSS %M kB" env OMP_NUM_THREADS=$NT $EXE < J$J.inp > J$J.out 2>&1 || true
  tail -1 J$J.out
  python3 $A/trove_energies.py J$J.out | tail -4
  if ! awk '/End of TROVE/{f=1} END{exit !f}' J$J.out; then echo "J=$J did not finish; stopping"; break; fi
done
