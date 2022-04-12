#!/bin/bash

source activate Jesse_QM_MM

export PATH=/home/mcdanielgroup/build/psi4/psi4/bin/bin:$PATH
export PSI_SCRATCH=/home/mcdanielgroup/data/scratch/
export PYTHONPATH=/home/mcdanielgroup/build/psi4/psi4/bin/lib:$PYTHONPATH

python run_QM_MM.py --pdb spce_spherical_cubic.pdb > output1
python run_QM_MM.py --pdb spce_spherical_orthorhombic.pdb > output2
