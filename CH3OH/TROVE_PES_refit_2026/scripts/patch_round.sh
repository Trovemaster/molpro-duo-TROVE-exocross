#!/bin/bash
# patch_round.sh <round> <prev_trove> : one hole-patching round
#   1. hole search (shell V_sep<30000) on the previous PEF
#   2. pseudo-points = random shell points + all holes found so far
#   3. Fortran refit (CH3OH_fitpes_ob.x)
set -e
R=$1; PREV=$2
ECAP=${ECAP:-30000}; TAG=${TAG:-W}
A=/home/user/ExoMol/CH3OH/pes_refit/analysis
F=/home/user/ExoMol/CH3OH/pes_refit/fit
cd $A
python3 -u shell_holes.py $PREV --ecap $ECAP --out holes${TAG}_round${R}.npz --seed $((7 + R)) > shell_holes${TAG}_round${R}.log 2>&1
tail -4 shell_holes${TAG}_round${R}.log
HOLES=$(ls holes${TAG}_round*.npz ${EXTRA_HOLES} 2>/dev/null | tr '\n' ' ')
python3 -u make_pseudo.py $PREV --holes $HOLES --ecap $ECAP --ratio 0.9 --n 1000 --seed $((11 + R)) \
  --out pseudo${TAG}_round$((R+1)).dat > make_pseudo${TAG}_round$((R+1)).log 2>&1
tail -3 make_pseudo${TAG}_round$((R+1)).log
mkdir -p $F/fit${TAG}_r$((R+1)) && cd $F/fit${TAG}_r$((R+1))
python3 $A/build_fit_input.py --terms $A/terms_E.txt --init $F/fitE_robust/00.res --itmax 8 \
  --nonlin-ivar 0 --term-ivar 1 --robust 1 --tol -1 --gamma 0.0005 --ps-emax 15000 --ps-in-python --ps-floor 0.003 \
  --ridge 1e-7 --pseudo $A/pseudo${TAG}_round$((R+1)).dat --trove-out CH3OH_PES_fit${TAG}_r$((R+1)).trove \
  --title "CH3OH PEF round $((R+1)): 4839 terms, PS+Watson, ridge 1e-7, pseudo-points" \
  --out CH3OH_fit${TAG}_r$((R+1)).inp > /dev/null
OMP_NUM_THREADS=48 $F/CH3OH_fitpes_ob.x < CH3OH_fit${TAG}_r$((R+1)).inp > CH3OH_fit${TAG}_r$((R+1)).out
grep -A10 "Unweighted rms" CH3OH_fit${TAG}_r$((R+1)).out
