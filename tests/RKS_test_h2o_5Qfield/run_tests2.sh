#!/bin/bash

module load intel/19.0
module load anaconda3

source activate gcc72

export PATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/bin:$PATH
export PSI_SCRATCH=~/scratch/
export PYTHONPATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/lib:$PYTHONPATH

python run_QM_MM.py --pme_grid 375 --cutoff 1.2 --quad_angle 2702 --quad_radial 89



