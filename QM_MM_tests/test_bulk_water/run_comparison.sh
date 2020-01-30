#!/bin/bash

module load intel/19.0
module load anaconda3

source activate gcc72

export PATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/bin:$PATH
export PSI_SCRATCH=~/scratch/
export PYTHONPATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/lib:$PYTHONPATH

# compare outputs :  run_QM_MM.py treats all MM water numerically with vext computed from PME, put into DFT quadrature
#                    run_QM_MM.py does the same calculation, but with three of the closest MM water molecules treated with analytic embedding.
#  this calculations should agree in the limit of a VERY LARGE PME grid, we find they agree well with a PME grid size of ~425
#  for a smaller PME grid, there will be too much error in the interpolation, and close MM molecules should be treated with analytic embedding
#
#  run_QM_MM3.py is same as run_QM_MM.py, except uses a different interpolation scheme (based on projecting to PME grid)

python run_QM_MM.py > output1
python run_QM_MM2.py > output2
python run_QM_MM3.py > output3
