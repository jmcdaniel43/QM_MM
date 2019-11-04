from __future__ import print_function
#*********** Psi4 Drivers
import psi4.core as core
#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
#*********** QM/MM classes
from QM_MM_classes import *
#**********************
# other stuff
from sys import stdout
from time import gmtime, strftime
from datetime import datetime



# *********************************************************************
#                     Create MM system object
#**********************************************************************

# Initialize: Input list of pdb and xml files
MMsys=MM( [ 'SC_start2_N2.pdb', ] , [ 'sapt_residues.xml' , 'graph_residue.xml' ] , [ 'sapt.xml', 'graph.xml' ] , [ 'sapt_Efield.xml' , 'graph.xml' ]   )
MMsys.set_periodic_residue(True)
# Initialze OpenMM API's
MMsys.set_platform('CPU')   # either 'CPU', 'OpenCL', or 'CUDA'
# Define MM constraint for N2
MMsys.setumbrella( 'N2', 'grph', 'C100', 2000.0 , 0.4 )   # molecule1, molecule2, atom2,  k (kJ/mol/nm^2) , r0 nm


# *********************************************************************
#                     Create QM system object
#**********************************************************************

# Define QM region and Initialize QM class
QMmolname='N2'
QMsys = QM( QMmolname )
# Fill QM region with atoms of QMmolname
QMsys.set_QM_region( MMsys )
# reference QM energy gas-phase
( QMsys.refenergy , QMsys.wfn ) = psi4.energy('scf/3-21G', return_wfn=True )


#**********************************************************************
#                     QM/MM Simulation
#**********************************************************************
# Get external field and positions from MMsys of QMregion
( positions, Efield ) = MMsys.get_positions_field( QMsys )

      
# update positions and external field on QMsys
QMsys.set_geometry( positions , 0 )
QMsys.set_external_field( *Efield ) 

# QM calculation
( QMsys.energy , QMsys.wfn ) = psi4.energy('scf/3-21G', return_wfn=True )

# QM charges
QMsys.get_DMA_charges()

# update MM with QM charges
MMsys.update_charges( QMsys )

# 10  MM integration steps
MMsys.simmd.step(10)


exit()
