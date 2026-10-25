#!/bin/bash
# final_checks.sh <final.trove> <final 00.res> <tag>
set -e
T=$(readlink -f $1); RES=$(readlink -f $2); TAG=$3
cd /home/user/ExoMol/CH3OH/pes_refit/analysis
echo "=== equilibrium";  python3 -u find_minimum.py $T | tail -13
echo "=== 1D/2D cuts + random box check"; python3 -u hole_check.py $T --no-min --nrand 200000
echo "=== shell holes V_sep<30000"; python3 -u shell_holes.py $T --ecap 30000 --nrand 1000000 --nmin 6000 --steps 200 --seed 303 --out holes_${TAG}_final30.npz | tail -5
echo "=== shell holes V_sep<40000"; python3 -u shell_holes.py $T --ecap 40000 --nrand 1000000 --nmin 6000 --steps 200 --seed 404 --out holes_${TAG}_final40.npz | tail -5
python3 plot_cuts.py $T cuts_1d_${TAG}.png
python3 plot_residuals.py $RES residuals_${TAG}.png
