from __future__ import print_function

import sys
sys.path.append('/home/mcdanielgroup/data/Jesse/QM_MM/lib/')
# append path to MM class library
sys.path.append('/home/mcdanielgroup/data/Jesse/Fixed_Voltage_OpenMM/lib/')
#********* import QMclass
from QM_classes import *
#********* import MMclass
from MM_classes import *

from routines import *
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
QMatoms_list=(0,)
# charge and spin
QMcharge=0
QMspin=1
#QMspin=2

# QM atoms and MM atoms with analytic Coulomb embedding
# note that QMregion should be ordered ( QMatoms , ... ), otherwise method set_geometry won't work properly.
# we currently don't have a check for this, as we anticipate that we will eventually generate QMregion automatically
# from QMatoms using a cutoff distance
QMregion_list=(0,)

# Make sure QMatoms is subset of QMregion
if not set(QMatoms_list).issubset(QMregion_list) :
   print(' QMatoms must be subset of QMregion !!')
   sys.exit()

# QMother is the difference between lists ..
QMother_list=np.setdiff1d( np.array( QMregion_list ) , np.array( QMatoms_list ) )
#**********************************************************************


DFT_functional='PBE'

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# Initialize: Input list of pdb and xml files, and QMregion
MMsys=MM( pdb_list = [ '../input_files/He_5ions.pdb', ] , residue_xml_list = [ '../input_files/sapt_residues.xml' , ] , ff_xml_list = [ '../input_files/sapt.xml', ] , QMregion_list = QMregion_list , cutoff = cutoff  )

# if periodic residue, call this
#MMsys.set_periodic_residue(True)

# set PME parameters.  This is important for control of accuracy for vext interpolation to DFT quadrature
# choice of alpha:  For n=43 grid, 60 Angstrom box, OpenMM chooses alpha= 2.389328 nm^-1
MMsys.setPMEParameters( pme_alpha=2.0 , pme_grid_a=pme_grid_size , pme_grid_b=pme_grid_size , pme_grid_c=pme_grid_size ) # alpha, nx, ny, nz

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform('Reference')   # only 'Reference' platform is currently implemented!

# IMPORTANT: generate exclusions for SAPT-FF
#sapt_exclusions = sapt_generate_exclusions(MMsys.simmd, MMsys.system, MMsys.modeller.positions)

# Umbrella potential on QM atoms
#MMsys.setumbrella( 'N2', 'grph', 'C100', 2000.0 , 0.4 )   # molecule1, molecule2, atom2,  k (kJ/mol/nm^2) , r0 nm


# *********************************************************************
#                     Create QM system object
#**********************************************************************

# Define QM region and Initialize QM class
# possible quadrature grids: see Psi4 manual

QMsys = QM( QMname = 'test' , basis = 'aug-cc-pvdz' , dft_spherical_points = quadrature_grid[0] , dft_radial_points = quadrature_grid[1] , scf_type = 'df' , qmmm='true' )

# get elements/charges of QM region atoms from MMsys ...
element_lists , charge_lists = MMsys.get_element_charge_for_atom_lists( [ QMatoms_list , QMother_list ] )
# Fill QM region with atoms.
QMsys.set_QM_region( element_lists , charge_lists , QMatoms_list, QMother_list )

# Get QM positions from MMsystem and set them in QMsys object
positions_lists = MMsys.get_positions_for_atom_lists([ QMatoms_list , QMother_list ] )
QMsys.set_QM_positions( positions_lists )

# set geometry of QM region
QMsys.set_geometry( charge = QMcharge, spin = QMspin )

#***************** test Psi4 quadrature **********************
#  input MMenv to Psi4, Psi4 will calculate the Coulomb potential
#  on the quadrature grid from this MMenv, and will use this
#  in the DFT quadrature.  This is for internal testing, so no field input from OpenMM...
#*************************************************************
MMenv = create_MM_env_full(  MMsys , QMatoms_list )

#print( 'MMenv ')
#print( MMenv )
#sys.exit()

# QM quadrature test---- note, if MMsys in input to Psi4 as **kwarg, then Psi4 internally evaluates Coulomb potential from MMenv on quadrature grid...
( QMsys.energy , QMsys.wfn ) = psi4.energy( DFT_functional , return_wfn=True , mm_env=MMenv )   # note, Psi4 will lowercase kwargs, use lowercase argument

print( 'QM energy ' , QMsys.energy )


exit()
