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

# file directories
pdb_dir = 'BMIM_box.pdb' 
res_dir = 'sapt_residues.xml'
ff_dir = 'sapt.xml'

# string arguments for MM and QM classes
platform = 'Reference'
functional = 'PBE'

# qmmm method arguments
qmmm_ewald = 1
qmmm_cutoffs = range( 3 , 37 , 3 )

# pme method rguments
pme_grid_size = 250
pme_alpha = 5.0

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit(2000)

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
QMatoms_list=(9900, 9901, 9902)

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# Initialize: Input list of pdb and xml files, and QMregion_list
MMsys=MM( pdb_list = [ pdb_dir , ] , residue_xml_list = [ res_dir , ] , ff_xml_list = [ ff_dir, ] , QMatoms_list = QMatoms_list , qmmm_ewald = qmmm_ewald )

# if periodic residue, call this
MMsys.set_periodic_residue( True )

# set PME parameters.  This is important for control of accuracy for vext interpolation to DFT quadrature
# choice of alpha:  For n=43 grid, 60 Angstrom box, OpenMM chooses alpha= 2.389328 nm^-1
MMsys.set_PMEParameters( pme_alpha=pme_alpha , pme_grid_a=pme_grid_size , pme_grid_b=pme_grid_size , pme_grid_c=pme_grid_size ) # alpha, nx, ny, nz

#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform( platform )   # only 'Reference' platform is currently implemented!

# *********************************************************************
#                     Create QM system object
#**********************************************************************

# Define QM region and Initialize QM class
# possible quadrature grids: see Psi4 manual
quadrature_grid = ( 302 , 75 )  # spherical points, radial points

QMsys = QM( QMname = 'test' , basis = 'aug-cc-pvdz' , dft_spherical_points = quadrature_grid[0] , dft_radial_points = quadrature_grid[1] , scf_type = 'df' , qmmm_ewald = qmmm_ewald , pme_grid_size = pme_grid_size , pme_alpha = pme_alpha )

#**********************************************************************
#                     QM/MM Simulation
#**********************************************************************

# instantiating list of energies to correspond with list of qmmm_cutoff distances
energies = []

print( 'calculating energy, force, and vext...' )

for i in range(len(qmmm_cutoffs)):

     # getting QMregion list
     MMsys.set_QMregion_parameters( QMatoms_list , qmmm_cutoffs[i] )
     MMsys.get_QMregion_list()

     # Make sure QMatoms_list is subset of QMregion_list
     if not set( QMatoms_list ).issubset( MMsys.QMregion_list ) :
          print( ' QMatoms_list must be subset of QMregion_list !!' )

     # QMother is the difference between lists ..
     QMother_list = np.setdiff1d( np.array( MMsys.QMregion_list ) , np.array( QMatoms_list ) )
     QMdrude_list = MMsys.QMdrudes_list

     # get elements/charges of QM region atoms from MMsys ...
     element_lists , charge_lists = MMsys.get_element_charge_for_atom_lists( [ QMatoms_list , QMother_list , QMdrude_list ] )

     # Fill QM region with atoms.
     QMsys.set_QM_region( element_lists , charge_lists , QMatoms_list, QMother_list , QMdrude_list )

     #******************* External potential on PME grid ******************
     if qmmm_ewald:

          # collecting state information
          state = MMsys.simmd.context.getState( getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True )
          print( 'done calculating energy, force, and vext' )

          # external potential on PME grid
          vext = state.getVext_grid()

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
          QMsys.get_PMEexclude( functional )

          # getting PME grid corrections and adding to vext_tot
          vext = QMsys.set_PMEexclude( real_pos , vext )

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
     if qmmm_ewald:
         ( QMsys.energy , QMsys.wfn ) = psi4.energy( functional , return_wfn=True , pme_grid_size=pme_grid_size , vexternal_grid=vext , box=box , interpolation_method="interpn" )
     else:
         ( QMsys.energy , QMsys.wfn ) = psi4.energy( functional , return_wfn=True )

     print( 'QM energy ' , QMsys.energy )

     energies.append( QMsys.energy )

# writing results to output file
fh = open( 'output' , 'w')

for d , e in zip( qmmm_cutoffs , energies ):
     fh.write( str(d) + ' ' + str(e) + '\n' )

fh.close()

sys.exit()
