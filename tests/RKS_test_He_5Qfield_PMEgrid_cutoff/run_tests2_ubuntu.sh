#!/bin/bash

source activate Jesse_QM_MM

export PATH=/home/mcdanielgroup/build/psi4/psi4/bin/bin:$PATH
export PSI_SCRATCH=/home/mcdanielgroup/data/scratch/
export PYTHONPATH=/home/mcdanielgroup/build/psi4/psi4/bin/lib:$PYTHONPATH

#job name
jobname="He_5Q"
# benchmark directory
benchmark_dir="reference_benchmark/"

# pairs of angular/radial quadrature settings
quad_angular=( 2702  )
quad_radial=(    89 )

# pairs of PMEgrid/cutoff settings
PME_grid=( 65  85 )
cutoff=( 1.2 0.8 ) # in nanometers ...

i=0
while [ $i  -lt ${#quad_angular[@]} ] ; do
   # find reference calculation
   file="${benchmark_dir}${jobname}_${quad_angular[$i]}_${quad_radial[$i]}.out"

   if [ -f ${file} ] ; then
      #echo ${file}
      # get HOMO, LUMO, and SCF reference energy
      homo_reference=`grep -A2 "Doubly Occupied" ${file} | grep "A" | cut -c 10-`
      lumo_reference=`grep -A2 "Virtual" ${file} | grep "A" | cut -c 10-25`
      scf_reference=`grep "Total Energy =" ${file} | cut -c 40-`
   else
      echo "can't find reference calculation ${file}"
   fi

   # loop over PMEgrid and cutoff settings
   j=0
   while [ $j  -lt ${#PME_grid[@]} ] ; do
       echo "test for PMEgrid : ${PME_grid[$j]} , cutoff : ${cutoff[$j]} , quadrature angular : ${quad_angular[$i]} radial:  ${quad_radial[$i]} ... " ;
       # now run test calculation
       python run_QM_MM.py --pme_grid ${PME_grid[$j]} --cutoff ${cutoff[$j]} --quad_angle ${quad_angular[$i]} --quad_radial ${quad_radial[$i]}  > testout

       # get HOMO, LUMO, and SCF test energy
       homo_test=`grep -A2 "Doubly Occupied" testout | grep "A" | cut -c 10-`
       lumo_test=`grep -A2 "Virtual" testout | grep "A" | cut -c 10-25`
       scf_test=`grep "Total Energy =" testout | cut -c 40-`
   
       echo " HOMO(test) HOMO(ref) LUMO(test) LUMO(ref) SCFenergy(test) SCFenergy(ref) "
       echo " ${homo_test} ${homo_reference} ${lumo_test} ${lumo_reference} ${scf_test} ${scf_reference} "

       ((j+=1));
   done
   ((i+=1)) ;
done


