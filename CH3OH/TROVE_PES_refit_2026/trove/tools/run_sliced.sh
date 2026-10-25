#!/bin/bash
# run_sliced.sh <rundir> <erange> <jmax> [threads] [first_step]
# Sliced TROVE J=0-contraction pipeline (see make_slice_inputs.py):
#   J0s -> (slice_s, conv_s, rm matelem s) for s=1..12 -> J1c..J<jmax>c
set -u
RUN=$1; ER=$2; JMAX=$3; NT=${4:-12}; FIRST=${5:-0}
B=/home/user/ExoMol/CH3OH/pes_refit/trove/bends
EXE=/home/user/ExoMol/CH3OH/pes_refit/trove/src_ob/j-trove-ob4.x
A=/home/user/ExoMol/CH3OH/pes_refit/analysis
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT OMP_STACKSIZE=1G
ulimit -s unlimited
cd $RUN
ok() { awk '/End of TROVE|\.\.\.done!/{f=1} /Program received signal|STOP|ERROR/{e=1} END{exit !(f && !e)}' $1; }
step() { local inp=$1; local out=${1%.inp}.out
  /usr/bin/time -f "$(basename $inp): %e s, maxRSS %M kB" $EXE < $inp > $out 2>&1
  tail -1 $out; if ! ok $out; then echo "FAILED: $out"; tail -20 $out; exit 1; fi; df -h . | tail -1 | awk '{print "   disk free:", $4}'; }
if [ $FIRST -le 0 ]; then
  python3 $B/make_slice_inputs.py . $ER 0 $JMAX > /dev/null
  step J0s.inp
  python3 $A/trove_energies.py J0s.out | tail -4
fi
ZPE=$(awk '/Zero-point-energy is/{z=$NF} END{print z}' J0s.out)
python3 $B/make_slice_inputs.py . $ER $ZPE $JMAX
for s in $(seq 1 12); do
  [ $s -lt $FIRST ] && continue
  step slice_$s.inp
  ls -la matelem$s.chk | awk '{print "   matelem'$s'.chk", $5/1e9, "GB"}'
  step conv_$s.inp
  rm -f matelem$s.chk
done
for J in $(seq 1 $JMAX); do
  step J${J}c.inp
  python3 $A/trove_energies.py J${J}c.out | tail -4
done
