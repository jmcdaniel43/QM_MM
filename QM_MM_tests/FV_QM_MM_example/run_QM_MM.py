from __future__ import print_function
import sys
# append path to QM class library
sys.path.append('/home/mcdanielgroup/data/Jesse/QM_MM/lib/')
# append path to MM class library
sys.path.append('/home/mcdanielgroup/data/Jesse/Fixed_Voltage_OpenMM/lib/')
#********* import QMclass
from QM_classes import *
#********* import MMclass
from MM_classes import *
from Fixed_Voltage_routines import *
from add_customnonbond_xml import add_CustomNonbondedForce_SAPTFF_parameters
#***************************
import numpy as np
import urllib.request
# other stuff
import sys
from sys import stdout
from time import gmtime, strftime
from datetime import datetime

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit(2000)

# *********************************************************************
#                     Define QM region
#
#  1) QMatoms_list: These are atoms treated quantum mechanically 
#
#  2) QMregion_list: In addition to quantum mechanical atoms, this includes MM atoms that should be treated
#     with analytic Coulomb embedding in Psi4.  In general QMatoms_list should be a subset of QMregion_list.
#
#  3) MMatoms are everything else.  Any atoms not in QMregion_list will contribute to vext_tot computed by OpenMM
#
#  4) In, general, should use at least a 3-5 Angstrom buffer of atoms around QMatoms_list that is included in QMregion_list.
#     these MM atoms will have analytic electrostatic embedding within Psi4.  The reason this is important is 2-fold.
#     First, even a small DFT quadrature grid will extend 3-4 Angstrom away from the nucleus.  This will overlap with
#     some atoms in the MM region, and so Coulomb interactions with these grid points will diverge.  Obviously this is not
#     a problem for analytic integrals, but is a problem for the numerical vext_tot, which may give 'infty' from Coulomb divergence
#     in OpenMM.
#     Second, this should significantly help errors due to interpolation from the PME_grid to the quadrature grid.  Without this buffer region,
#     interpolation errors will be significant, especially for nuclear energy contribution
#  
#
hartree_to_kjmol = 2625.4996
nm_to_bohr = 18.89726

# these are atom indices from input '*.pdb' file that define the atoms in the QM region
# remember, indexing starts at zero ...

# QM atoms only
QMatoms_list=(8152,8153,8154,8155,8156,8157)
# charge and spin
QMcharge=0
QMspin=1
#QMspin=2

# QM atoms and MM atoms with analytic Coulomb embedding
QMregion_list=(8152,8153,8154,8155,8156,8157)

# Make sure QMatoms_list is subset of QMregion_list
if not set(QMatoms_list).issubset(QMregion_list) :
   print(' QMatoms_list must be subset of QMregion_list !!')
   sys.exit()

# QMother is the difference between lists ..
QMother_list=np.setdiff1d( np.array( QMregion_list ) , np.array( QMatoms_list ) )
#**********************************************************************

#************************** download SAPT-FF force field files from github
url1 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt.xml"
url2 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt_residues.xml"
#url3 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt_exclusions.py"
filename1, headers = urllib.request.urlretrieve(url1, "ffdir/sapt.xml")
filename2, headers = urllib.request.urlretrieve(url2, "ffdir/sapt_residues.xml")
#filename3, headers = urllib.request.urlretrieve(url3, "sapt_exclusions.py")

# add extra CustomNonbonded force parameters to .xml file
add_CustomNonbondedForce_SAPTFF_parameters( xml_base = "ffdir/sapt.xml" , xml_param = "ffdir/graph_customnonbonded.xml" , xml_combine = "ffdir/sapt_add.xml" )


DFT_functional='PBE'
pme_grid_size=200

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# set applied voltage in Volts
Voltage = 4.0 # in Volts, units will be internally converted later...


# electrode names used to exclude intra-electrode non-bonded interactions ...
#cathode_name="cath"; anode_name = "anod"
# use chain indices instead of residue names to identify electrode...
cathode_index=(0,2); anode_index=(1,3) # note chain indices start at 0 ...

# Initialize: Input list of pdb and xml files, and QMregion_list

MMsys=MM( pdb_list = [ 'graphene_BMIM_BF4_acnt.pdb', ] , residue_xml_list = [ './ffdir/sapt_residues.xml' , './ffdir/graph_residue_c.xml', './ffdir/graph_residue_n.xml' ] , ff_xml_list = [ './ffdir/sapt_add.xml', './ffdir/graph_c_freeze.xml', './ffdir/graph_n_freeze.xml' , './ffdir/graph.xml' ] , QMregion_list = QMregion_list , cutoff = 1.0*nanometer )

# if periodic residue, call this
MMsys.set_periodic_residue(True)

# set PME parameters.  This is important for control of accuracy for vext interpolation to DFT quadrature
# choice of alpha:  For n=43 grid, 60 Angstrom box, OpenMM chooses alpha= 2.389328 nm^-1
MMsys.setPMEParameters( pme_alpha=2.4 , pme_grid_a=pme_grid_size , pme_grid_b=pme_grid_size , pme_grid_c=pme_grid_size ) # alpha, nx, ny, nz

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform('Reference')   # only 'Reference' platform is currently implemented!

# initialize Virtual Electrodes, these are electrode `sheets' that solve electrostatics for constant Voltage ...
# can choose electrodes by residue name (this is default)
# can also choose electrodes by chain name (set chain=True)
# can input tuple exclude_element with elements to exclude from virtual electrode, such as dummy Hydrogen atoms ...
#MMsys.initialize_electrodes( Voltage, cathode_identifier = cathode_name , anode_identifier = anode_name , chain=False , exclude_element=("H",) )  # use residue names as identifiers ...
MMsys.initialize_electrodes( Voltage, cathode_identifier = cathode_index , anode_identifier = anode_index , chain=True , exclude_element=("H",) )  # chain indices instead of residue names

# initialize atoms indices of electrolyte, we need this for analytic charge correction.  Currently we electrode residue > 100 atoms, electrolyte residue < 100 atoms ... this should be fine?
MMsys.initialize_electrolyte(Natom_cutoff=100)  # make sure all electrode residues have greater than, and all electrolyte residues have less than this number of atoms...


# IMPORTANT: generate exclusions for SAPT-FF.  If flag_SAPT_FF_exclusions = True , then will assume SAPT-FF force field and put in appropriate exclusions.
# set flag_SAPT_FF_exclusions = False if not using SAPT-FF force field
MMsys.generate_exclusions( flag_SAPT_FF_exclusions = True )

# Umbrella potential on QM atoms
#MMsys.setumbrella( 'N2', 'grph', 'C100', 2000.0 , 0.4 )   # molecule1, molecule2, atom2,  k (kJ/mol/nm^2) , r0 nm

# Fixed Voltage Electrostatics ...
MMsys.Poisson_solver_fixed_voltage( Niterations=3 )


# *********************************************************************
#                     Create QM system object
#**********************************************************************

# Define QM region and Initialize QM class
# possible quadrature grids: see Psi4 manual
#quadrature_grid = ( 50 , 12 )  # spherical points, radial points
#quadrature_grid = ( 302 , 75 )  # spherical points, radial points
quadrature_grid = ( 2702 , 89 )  # spherical points, radial points

QMsys = QM( QMname = 'test' , basis = '6-311G' , dft_spherical_points = quadrature_grid[0] , dft_radial_points = quadrature_grid[1] , scf_type = 'df' , qmmm='true' )

# get elements/charges of QM region atoms from MMsys ...
element_lists , charge_lists = MMsys.get_element_charge_for_atom_lists( [ QMatoms_list , QMother_list ] )
# Fill QM region with atoms.
QMsys.set_QM_region( element_lists , charge_lists , QMatoms_list, QMother_list )

#**********************************************************************
#                     QM/MM Simulation
#**********************************************************************

print( 'calculating energy, force, and vext...')
#******************* External potential on PME grid ******************
#state = MMsys.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True)
state = MMsys.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True)
print('done calculating energy, force, and vext')

# external potential on PME grid
vext_tot = state.getVext_grid()
# PME grid positions
#PME_grid_positions = state.getPME_grid_positions()      

#** vext_tot from OpenMM is in kJ/mol/e  units, convert to Hartree for input to Psi4
vext_tot = np.array( vext_tot ) / hartree_to_kjmol

#** PME_grid_positions from OpenMM is in nanometers, convert to Bohr for input to Psi4
#PME_grid_positions = np.array( PME_grid_positions ) * nm_to_bohr

#** box vectors in units of Bohr
box = get_box_vectors_Bohr( state , nm_to_bohr )

# Get QM positions from MMsystem and set them in QMsys object
positions_lists = MMsys.get_positions_for_atom_lists([ QMatoms_list , QMother_list ] )
QMsys.set_QM_positions( positions_lists )

# update positions on QMsys
QMsys.set_geometry( charge = QMcharge , spin = QMspin )


# QM calculation
#*********************************************************
#       2 ways to intepolate vext from PME grid to DFT quadrature.  These are controlled by input **kwarg to psi4.energy
#  ****** Method 1:  Input PME_grid_positions in real space, interpolate in real space *********
#        LIMITATIONS- a) this doesn't consider PBC, if any points in quadrature grid fall outside of principle simulation box, then interpolation will crash!
#                     b) this can't handle non-cubic (triclinic) boxes.  Use Method 2 if either a) or b) is a problem    
# 2 options for scipy interpolation: set interpolation_method = "interpn" or "griddata".  "interpn" should be much faster and is for regularly spaced grids
#( QMsys.energy , QMsys.wfn ) = psi4.energy( DFT_functional , return_wfn=True , pme_grid_size=pme_grid_size , vexternal_grid=vext_tot , pmegrid_xyz = PME_grid_positions , interpolation_method="interpn" )
#
#  ******  Method 2:  Input box vectors, project DFT quadrature grid to PME grid, and interpolate on PME grid
#         this works with arbitrary triclinic boxes, and imposes PBC on DFT quadrature grid
( QMsys.energy , QMsys.wfn ) = psi4.energy( DFT_functional , return_wfn=True , pme_grid_size=pme_grid_size , vexternal_grid=vext_tot , box = box , interpolation_method="interpn" )

print( 'QM energy ' , QMsys.energy )


exit()
