from __future__ import print_function
#*********** Psi4 Drivers
import psi4.core as core
#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
#*********** QM/MM classes
from QM_MM_classes import *
#******** this is module that goes with sapt force field files to generate exclusions
from sapt_exclusions import *
#***************************
#from routines import trim_shift_PME_grid
import numpy as np
# other stuff
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pme_grid", type=str, help="pme_grid_size")
parser.add_argument("--cutoff", type=str, help="real space cutoff")
parser.add_argument("--quad_angle", type=str, help="angular quadrature grid")
parser.add_argument("--quad_radial", type=str, help="radial quadrature grid")
args = parser.parse_args()

pme_grid_size=int(args.pme_grid)
cutoff = float(args.cutoff) * nanometer
quadrature_grid = ( int(args.quad_angle) , int(args.quad_radial) )

# *********************************************************************
#                     Define QM region
#
#  1) QMatoms: These are atoms treated quantum mechanically 
#
#  2) QMregion: In addition to quantum mechanical atoms, this includes MM atoms that should be treated
#     with analytic Coulomb embedding in Psi4.  In general QMatoms should be a subset of QMregion.
#
#  3) MMatoms are everything else.  Any atoms not in QMregion will contribute to vext_tot computed by OpenMM
#
#  4) In, general, should use at least a 3-5 Angstrom buffer of atoms around QMatoms that is included in QMregion.
#     these MM atoms will have analytic electrostatic embedding within Psi4.  The reason this is important is 2-fold.
#     First, even a small DFT quadrature grid will extend 3-4 Angstrom away from the nucleus.  This will overlap with
#     some atoms in the MM region, and so Coulomb interactions with these grid points will diverge.  Obviously this is not
#     a problem for analytic integrals, but is a problem for the numerical vext_tot, which may give 'infty' from Coulomb divergence
#     in OpenMM.
#     Second, this should significantly help errors due to interpolation from the PME_grid to the quadrature grid.  Without this buffer region,
#     interpolation errors will be significant, especially for nuclear energy contribution
#  
#

# these are atom indices from input '*.pdb' file that define the atoms in the QM region
# remember, indexing starts at zero ...
#QMatoms=(2,3)

# QM atoms only
QMatoms=(0,)
# charge and spin
QMcharge=0
QMspin=1
#QMspin=2

# QM atoms and MM atoms with analytic Coulomb embedding
# note that QMregion should be ordered ( QMatoms , ... ), otherwise method set_geometry won't work properly.
# we currently don't have a check for this, as we anticipate that we will eventually generate QMregion automatically
# from QMatoms using a cutoff distance
QMregion=(0,1,2,3)

# Make sure QMatoms is subset of QMregion
if not set(QMatoms).issubset(QMregion) :
   print(' QMatoms must be subset of QMregion !!')
   sys.exit()

#**********************************************************************


DFT_functional='PBE'

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# Initialize: Input list of pdb and xml files, and QMregion
MMsys=MM( pdb_list = [ '../input_files/He_5ions.pdb', ] , residue_xml_list = [ '../input_files/sapt_residues.xml' , ] , ff_xml_list = [ '../input_files/sapt.xml', ] , QMregion = QMregion , cutoff = cutoff  )

# if periodic residue, call this
#MMsys.set_periodic_residue(True)

# set PME parameters.  This is important for control of accuracy for vext interpolation to DFT quadrature
# choice of alpha:  For n=43 grid, 60 Angstrom box, OpenMM chooses alpha= 2.389328 nm^-1
MMsys.setPMEParameters( pme_alpha=2.0 , pme_grid_a=pme_grid_size , pme_grid_b=pme_grid_size , pme_grid_c=pme_grid_size ) # alpha, nx, ny, nz

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform('Reference')   # only 'Reference' platform is currently implemented!

# IMPORTANT: generate exclusions for SAPT-FF
sapt_exclusions = sapt_generate_exclusions(MMsys.simmd, MMsys.system, MMsys.modeller.positions)

# Umbrella potential on QM atoms
#MMsys.setumbrella( 'N2', 'grph', 'C100', 2000.0 , 0.4 )   # molecule1, molecule2, atom2,  k (kJ/mol/nm^2) , r0 nm


# *********************************************************************
#                     Create QM system object
#**********************************************************************

# Define QM region and Initialize QM class
# possible quadrature grids: see Psi4 manual
#quadrature_grid = ( 50 , 12 )  # spherical points, radial points
#quadrature_grid = ( 302 , 75 )  # spherical points, radial points
#quadrature_grid = ( 2702 , 89 )  # spherical points, radial points

QMsys = QM( QMname = 'test' , basis = 'aug-cc-pvdz' , dft_spherical_points = quadrature_grid[0] , dft_radial_points = quadrature_grid[1] , scf_type = 'df' , qmmm='true' )

# Fill QM region with atoms.  Use MMsys to get element types
QMsys.set_QM_region( MMsys, QMregion, QMatoms )

# initial QM energy, Get QM positions from MMsystem
positions  = MMsys.get_positions_QM( QMsys )

# set geometry of QM region
QMsys.set_geometry( positions , charge = QMcharge, spin = QMspin )


#**********************************************************************
#                     QM/MM Simulation
#**********************************************************************
# Get QM positions from MMsystem
positions  = MMsys.get_positions_QM( QMsys )

#******************* External potential on PME grid ******************
state = MMsys.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True)

# external potential on PME grid
vext_tot = state.getVext_grid()
# PME grid positions
PME_grid_positions = state.getPME_grid_positions()      
#** vext_tot from OpenMM is in kJ/mol/e  units, convert to Hartree for input to Psi4
vext_tot = np.array( vext_tot ) / 2625.4996
#** PME_grid_positions from OpenMM is in nanometers, convert to Bohr for input to Psi4
PME_grid_positions = np.array( PME_grid_positions ) * 18.89726

# update positions on QMsys
QMsys.set_geometry( positions , charge = QMcharge , spin = QMspin )

# QM calculation
# 2 options for scipy interpolation: set interpolation_method = "interpn" or "griddata".  "interpn" should be much faster and is for regularly spaced grids
( QMsys.energy , QMsys.wfn ) = psi4.energy( DFT_functional , return_wfn=True , pme_grid_size=pme_grid_size , vexternal_grid=vext_tot , pmegrid_xyz = PME_grid_positions , interpolation_method="interpn" )

print( 'QM energy ' , QMsys.energy )

sys.exit()

# QM charges
QMsys.get_DMA_charges()

# update MM with QM charges
MMsys.update_charges( QMsys )

# 10  MM integration steps
MMsys.simmd.step(10)


exit()
