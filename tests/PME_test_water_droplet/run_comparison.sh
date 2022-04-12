#!/bin/bash

module load intel/19.0
module load anaconda3

source activate gcc72

export PATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/bin:$PATH
export PSI_SCRATCH=~/scratch/
export PYTHONPATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/lib:$PYTHONPATH

python run_test_vext.py > output
