#!/bin/bash
# chain_after_j0s.sh <rundir> <erange> <jmax> <threads>: wait for <rundir>/J0s.out to finish, then run
# the sliced J=0-contraction pipeline (slices 1..12, conversions, J=1..jmax).
RUN=$1; ER=$2; JMAX=$3; NT=$4
B=/home/user/ExoMol/CH3OH/pes_refit/trove/bends
while true; do
  if awk '/End of TROVE/{f=1} END{exit !f}' $RUN/J0s.out 2>/dev/null; then break; fi
  if awk '/Program received signal|STOP|Killed/{f=1} END{exit !f}' $RUN/J0s.out 2>/dev/null; then echo "J0s failed in $RUN"; exit 1; fi
  sleep 30
done
echo "J0s finished in $RUN at $(date)"
$B/run_sliced.sh $RUN $ER $JMAX $NT 1
