#!/bin/bash

module load intel/19.0
module load anaconda3

source activate gcc72

export PATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/bin:$PATH
export PSI_SCRATCH=~/scratch/
export PYTHONPATH=/nv/hp22/jmcdaniel43/data/Programs/psi4/psi4/bin/lib:$PYTHONPATH

#job name
jobname="He_5Q"
# benchmark directory
benchmark_dir="reference_benchmark/"

# pairs of angular/radial quadrature settings
quad_angular=( 50  302  2702  5810 )
quad_radial=( 12  75  89  131 )


i=0
while [ $i  -lt ${#quad_angular[@]} ] ; do
   echo "test for quadrature angular : ${quad_angular[$i]} radial:  ${quad_radial[$i]} ... " ;
   # find reference calculation
   file="${benchmark_dir}${jobname}_${quad_angular[$i]}_${quad_radial[$i]}.out"

   if [ -f ${file} ] ; then
      echo ${file}
      # get HOMO, LUMO, and SCF reference energy
      homo_reference=`grep -A2 "Doubly Occupied" ${file} | grep "A" | cut -c 10-`
      lumo_reference=`grep -A2 "Virtual" ${file} | grep "A" | cut -c 10-25`
      scf_reference=`grep "Total Energy =" ${file} | cut -c 40-`
   else
      echo "can't find reference calculation ${file}"
   fi

   # now run test calculation
   # note pme_grid and cutoff setting are not used here, as MMfield on quadrature grid is computed internally by Psi4 ...
   python run_QM_MM.py --pme_grid 60 --cutoff 1.2 --quad_angle ${quad_angular[$i]} --quad_radial ${quad_radial[$i]}  > testout

   # get HOMO, LUMO, and SCF test energy
   homo_test=`grep -A2 "Doubly Occupied" testout | grep "A" | cut -c 10-`
   lumo_test=`grep -A2 "Virtual" testout | grep "A" | cut -c 10-25`
   scf_test=`grep "Total Energy =" testout | cut -c 40-`
   
   echo " HOMO(test) HOMO(ref) LUMO(test) LUMO(ref) SCFenergy(test) SCFenergy(ref) "
   echo " ${homo_test} ${homo_reference} ${lumo_test} ${lumo_reference} ${scf_test} ${scf_reference} "

   ((i+=1)) ;
done


