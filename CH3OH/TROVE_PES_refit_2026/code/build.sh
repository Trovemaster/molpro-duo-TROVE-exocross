#!/bin/bash
# build.sh: compile CH3OH_fitpes (netlib LAPACK) and CH3OH_fitpes_ob (OpenBLAS ILP64)
set -e
cd "$(dirname "$0")"
OB=${OB:-$HOME/.local/lib/python3.13/site-packages/numpy.libs}
OBLIB=$(ls $OB/libscipy_openblas64_*.so 2>/dev/null | head -1)
gfortran -cpp -O3 -march=native -fopenmp CH3OH_fitpes.f90 -llapack -lblas -o CH3OH_fitpes.x
if [ -n "$OBLIB" ]; then
  gfortran -cpp -DOPENBLAS64 -O3 -march=native -fopenmp CH3OH_fitpes.f90 -o CH3OH_fitpes_ob.x \
    -L$OB -l:$(basename $OBLIB) -Wl,-rpath,$OB -llapack -lblas
fi
ls -la CH3OH_fitpes*.x
