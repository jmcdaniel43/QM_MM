from __future__ import print_function
import sys

# append QM class library to path
sys.path.append('/nv/hp22/jpederson6/data/fv_qmmm/QM_MM/lib/')

# append MM class library to path
sys.path.append('/nv/hp22/jpederson6/data/fv_qmmm/Fixed_Voltage_OpenMM/lib/')

#********* import QMclass
from QM_classes import *

#********* import MMclass
from MM_classes import *

# other stuff
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import numpy as np
import argparse

# collecting and parsing user arguments
parser = argparse.ArgumentParser( description='cutoff distance in angstroms and output file name' )

# file directories
parser.add_argument( 'pdb_dir' , type=str )
parser.add_argument( 'res_dir' , type=str )
parser.add_argument( 'ff_dir' , type=str )
parser.add_argument( 'out_dir' , type=str )

# string arguments for MM and QM classes
parser.add_argument( 'platform' , type=str )
parser.add_argument( 'functional' , type=str )

# qmmm method arguments
parser.add_argument( 'qmmm_ewald' , type=int )
parser.add_argument( 'qmmm_cutoff' , type=float )

# pme method rguments
parser.add_argument( 'pme_grid_size' , type=int )
parser.add_argument( 'pme_alpha' , type=float )

args = parser.parse_args()

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit(2000)

# DEPRECATED - defined in QMsys
# important constants
#hartree_to_kjmol = 2625.4996
#nm_to_bohr = 18.89726

# setting charge and spin
QMcharge=0
QMspin=1
#QMspin=2

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

# these are atom indices from input '*.pdb' file that define the atoms in the QM region
# remember, indexing starts at zero ...

# QM atoms only

#QMatoms_list=(0,1,2)
#QMatoms_list=(670, 671, 672)
QMatoms_list=(9900, 9901, 9902)
#QMatoms_list=(668,669,670)
#QMatoms_list=(1524,1525,1526,1599,1600,1601,525,526,527,1833,1834,1835,168,169,170)
#QMatoms_list=(1524,1525,1526,1599,1600,1601,525,526,527)
#QMatoms_list=(1599,1600,1601,1833,1834,1835)
#QMatoms_list=(525,526,527)

#***********************DEPRECATED*************************************
# QM atoms and MM atoms with analytic Coulomb embedding
# we anticipate that we will eventually generate QMregion_list automatically
# from QMatoms_list using a cutoff distance
#QMregion_list=(1599,1600,1601,1833,1834,1835)
#QMregion_list=(1524,1525,1526,1599,1600,1601,525,526,527,1833,1834,1835,168,169,170)
#QMregion_list=(669, 670, 671, 66, 67, 68, 1281, 1282, 1283, 1392, 1393, 1394)
#**********************************************************************

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# Initialize: Input list of pdb and xml files, and QMregion_list
#MMsys=MM( pdb_list = [ 'md_finished.pdb', ] , residue_xml_list = [ 'sapt_residues.xml' , ] , ff_xml_list = [ 'sapt.xml', ] , QMatoms_list = QMatoms_list , QMMM_ewald = QMMM_ewald , QM_cutoff = args.QM_cutoff )
MMsys=MM( pdb_list = [ args.pdb_dir , ] , residue_xml_list = [ args.res_dir , ] , ff_xml_list = [ args.ff_dir, ] , QMatoms_list = QMatoms_list , qmmm_ewald = args.qmmm_ewald , qmmm_cutoff = args.qmmm_cutoff )
#MMsys=MM( pdb_list = [ 'spce.pdb', ] , residue_xml_list = [ 'h2o_res_pole.xml' , ] , ff_xml_list = [ 'h2o_pole.xml', ] , QMatoms_list = QMatoms_list , QMMM_ewald = QMMM_ewald , QM_cutoff = args.QM_cutoff )
#MMsys=MM( pdb_list = [ 'spce.pdb', ] , residue_xml_list = [ 'residues.xml' , ] , ff_xml_list = [ 'spce.xml', ] , QMatoms_list = QMatoms_list , QMMM_ewald = QMMM_ewald , QM_cutoff = args.QM_cutoff )

# if periodic residue, call this
MMsys.set_periodic_residue( True )

# set PME parameters.  This is important for control of accuracy for vext interpolation to DFT quadrature
# choice of alpha:  For n=43 grid, 60 Angstrom box, OpenMM chooses alpha= 2.389328 nm^-1
MMsys.set_PMEParameters( pme_alpha=args.pme_alpha , pme_grid_a=args.pme_grid_size , pme_grid_b=args.pme_grid_size , pme_grid_c=args.pme_grid_size ) # alpha, nx, ny, nz

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform( args.platform )   # only 'Reference' platform is currently implemented!

# IMPORTANT: generate exclusions for SAPT-FF
#exclusions = electrode_sapt_generate_exclusions(MMsys.simmd, MMsys.system, MMsys.modeller.positions, cathode_name , anode_name)

# Umbrella potential on QM atoms
#MMsys.setumbrella( 'N2', 'grph', 'C100', 2000.0 , 0.4 )   # molecule1, molecule2, atom2,  k (kJ/mol/nm^2) , r0 nm

MMsys.get_QMregion_list()

# Make sure QMatoms_list is subset of QMregion_list
if not set( QMatoms_list ).issubset( MMsys.QMregion_list ) :
     print( ' QMatoms_list must be subset of QMregion_list !!' )

# QMother is the difference between lists ..
QMother_list = np.setdiff1d( np.array( MMsys.QMregion_list ) , np.array( QMatoms_list ) )
QMdrude_list = MMsys.QMdrudes_list

# *********************************************************************
#                     Create QM system object
#**********************************************************************

# Define QM region and Initialize QM class
# possible quadrature grids: see Psi4 manual
#quadrature_grid = ( 50 , 12 )  # spherical points, radial points
#quadrature_grid = ( 302 , 50 )  # spherical points, radial points
quadrature_grid = ( 302 , 75 )  # spherical points, radial points
#quadrature_grid = ( 2702 , 89 )  # spherical points, radial points

QMsys = QM( QMname = 'test' , basis = 'aug-cc-pvdz' , dft_spherical_points = quadrature_grid[0] , dft_radial_points = quadrature_grid[1] , scf_type = 'df' , qmmm_ewald = args.qmmm_ewald , pme_grid_size = args.pme_grid_size , pme_alpha = args.pme_alpha )

# get elements/charges of QM region atoms from MMsys ...
element_lists , charge_lists = MMsys.get_element_charge_for_atom_lists( [ QMatoms_list , QMother_list , QMdrude_list ] )

# Fill QM region with atoms.
QMsys.set_QM_region( element_lists , charge_lists , QMatoms_list, QMother_list , QMdrude_list )

#**********************************************************************
#                     QM/MM Simulation
#**********************************************************************

print( 'calculating energy, force, and vext...' )

#******************* External potential on PME grid ******************
if args.qmmm_ewald:

    # collecting state information
    state = MMsys.simmd.context.getState( getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True )
    print( 'done calculating energy, force, and vext' )

    # external potential on PME grid
    vext = state.getVext_grid()

    # PME grid positions
    #PME_grid_positions = state.getPME_grid_positions()      

    #** PME_grid_positions from OpenMM is in nanometers, convert to Bohr for input to Psi4
    #PME_grid_positions = np.array( PME_grid_positions ) * nm_to_bohr

    #** box vectors in units of Bohr
    box = get_box_vectors_Bohr( state )

    # Get QM positions from MMsystem and set them in QMsys object
    positions_lists , real_pos = MMsys.get_positions_for_atom_lists([ QMatoms_list , QMother_list , QMdrude_list] )
    QMsys.set_QM_positions( positions_lists )

    # set geometry of QM region
    QMsys.set_geometry( charge = QMcharge, spin = QMspin )

    # setting PME grid positions
    QMsys.set_PMEgrid_xyz( box )

    # getting PME grid points to exclude from reciprocal space
    QMsys.get_PMEexclude( args.functional )

    # getting PME grid corrections and adding to vext_tot
    vext = QMsys.set_PMEexclude( real_pos , vext )

    # kind of really hcky right now
    #all_pos = state.getPositions()
    #MM_pos = []
    #for i in all_pos:
    #    MM_pos.append(i._value)
    #MM_pos = np.array(MM_pos)
    #index = range(len(MM_pos))
    #QM_index = np.concatenate((np.array(QMatoms_list),np.array(QMother_list),np.array(QMdrude_list)))
    #MM_pos = np.delete( MM_pos , QM_index , 0 )
    #index = np.delete( np.array(index) , QM_index )
    #element_lists , charge_lists = MMsys.get_element_charge_for_atom_lists( [index.tolist()] )
    #QMsys.set_real_correction( args.pme_alpha , MM_pos , vext , charge_lists[0] )

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
if args.qmmm_ewald:
    ( QMsys.energy , QMsys.wfn ) = psi4.energy( args.functional , return_wfn=True , pme_grid_size=args.pme_grid_size , vexternal_grid=vext , box=box , interpolation_method="interpn" )
else:
    ( QMsys.energy , QMsys.wfn ) = psi4.energy( args.functional , return_wfn=True )
print( 'QM energy ' , QMsys.energy )

fh = open( args.out_dir + '.out' , 'w')
fh.write( str(QMsys.energy) )
fh.close()

#fh = open( args.output + '.out' , 'w')
#for i in range(len(vext1)):
#    fh.write( str(vext1[i]) + '\t\t' + str(vext[i]) + '\n' )
#fh.close()

sys.exit()
