from sys import stdout

import psi4
import psi4.core as core

from psi4.driver import *

from scipy.special import erf
from scipy.special import erfc

from psi4.driver.procrouting.qmmm.grid_interface import project_to_PME_grid
from psi4.driver.procrouting.qmmm.grid_interface import calculate_inverse_box_vectors

import sys

import time
import subprocess

# conversions
nm_to_bohr = 18.89726
hartree_to_kjmol = 2625.4996

#*************************** README  **************************************
#  This module defines QM classes that is used in a QM/MM simulation
#  Interfacing Psi4  and OpenMM
#
#  Because these codes use different units and datastructures, conversions
#  between these are done within the QM/MM classes, so that input and output
#  passed to/from the classes is given in units/datastructures of the
#  particular libraries (Psi4, OpenMM) that it comes from/goes to. 
#
#**************************************************************************

#  This class holds atom info for QMregion.  Not to be confused with atom class of OpenMM ...
class atom_QM(object):
    def __init__(self, element, charge, atom_index ):
        self.element = element
        self.charge  = charge
        self.atom_index = atom_index
        # initialize positions to zero in lieu of anything better...
        self.x = 0.0; self.y=0.0; self.z=0.0

    def set_xyz( self, x , y , z ):
        self.x = x
        self.y = y 
        self.z = z


#  This class is used to initialize a QM simulation
class QM(object):
    def __init__( self , **kwargs ):

          # setup options for DMA
          radii=['H',0.53,'C',0.53,'N',0.53,'O',0.53]

          core.set_local_option('GDMA', 'GDMA_RADIUS', radii )
          core.set_local_option('GDMA', 'GDMA_LIMIT', 4 )
          core.set_local_option('GDMA', 'GDMA_MULTIPOLE_UNITS', 'AU' )
          core.set_local_option('GDMA', 'GDMA_SWITCH', 0.0 )

          #storing inputs for later
          self.inputs = kwargs

          # reading inputs from **kwargs
          if 'out_dir' in kwargs :
              self.molname = kwargs['out_dir'].split()[0]
          if 'basis_set' in kwargs :
              self.basis = kwargs['basis_set']              
          if 'functional' in kwargs :
              self.DFT_functional = kwargs['functional']              
          if 'quadrature_spherical' in kwargs :
              self.dft_spherical_points = int(kwargs['quadrature_spherical'])
          if 'quadrature_radial' in kwargs :
              self.dft_radial_points = int(kwargs['quadrature_radial'])
          if 'qmmm_ewald' in kwargs :
              self.qmmm_ewald = eval(kwargs['qmmm_ewald'])
          if 'qmmm_tare' in kwargs :
              self.qmmm_tare = eval(kwargs['qmmm_tare'])
          if 'QMcharge' in kwargs :
              self.QMcharge = int(kwargs['QMcharge'])
          if 'QMspin' in kwargs :
              self.QMspin = int(kwargs['QMspin'])
          if 'pme_real_correction' in kwargs :
              self.pme_real_correction = eval(kwargs['pme_real_correction'])

          # setting scf_type
          self.scf_type = 'df'

          # set basis , quadrature settings, scf_type, qmmm ...
          psi4.set_options( {'basis' : self.basis , 'dft_spherical_points' : self.dft_spherical_points , 'dft_radial_points' : self.dft_radial_points , 'scf_type' : self.scf_type , 'qmmm' : str(self.qmmm_ewald).lower() } )

          # setting pme options
          if self.qmmm_ewald:
              self.pme_grid_size = int(kwargs['pme_grid_size'])
              self.pme_alpha = float(kwargs['pme_alpha'])

    #**********************************
    # Reinitializes the QM system object
    # kwargs is a dictionary of arguments
    #**********************************
    def reset( self , **kwargs ):
          self.__init__( **kwargs )

    #**********************************
    # QMatoms_list is tuple specifying atom indices of QM atom
    # QMother_list specifies additional MM embedding atoms
    # element_lists / charge_lists correspond to [ QMatoms_list , QMother_list ] ...
    #**********************************
    def set_QM_region(self, element_lists , charge_lists , QMatoms_list , QMother_list , QMdrude_list ):
          # create lists of atom_QM objects: there are two lists here:
          # QMatoms stores data for QMatoms
          # QMother stores data for MM embedding atoms
          self.QMatoms = []
          self.QMother = []
          self.QMdrude = []

          # first create self.QMatoms ..
          for i_atom in range(len(QMatoms_list)):
              # create atom_QM object
              atom_i = atom_QM( element_lists[0][i_atom] , charge_lists[0][i_atom] , QMatoms_list[i_atom] )
              self.QMatoms.append( atom_i )  #QMatom

         # now create self.QMother ..
          for i_atom in range(len(QMother_list)):
              # create atom_QM object
              atom_i = atom_QM( element_lists[1][i_atom] , charge_lists[1][i_atom] , QMother_list[i_atom] )
              self.QMother.append( atom_i )  #MM embedding atom

         # finally create self.QMdrude ..
          for i_atom in range(len(QMdrude_list)):
              # create atom_QM object
              atom_i = atom_QM( element_lists[2][i_atom] , charge_lists[2][i_atom] , QMdrude_list[i_atom] )
              self.QMdrude.append( atom_i )  #MM embedding atom


    #***************************************
    # this method sets positions in self.QMatom and self.QMother atom objects lists ...
    #
    # input positions_lists correspond to [ QMatoms_list , QMother_list ] ...
    #
    #   ***** Unit conversion *********
    #   OpenMM gives positions in nanometer
    #   Psi4 input positions in Angstrom
    #*****************************************
    def set_QM_positions( self, positions_lists ):
          lengthconv = 10.0  # nm to angstrom

          # first set positions in self.QMatoms ..
          for i_atom in range(len(self.QMatoms)):
              atom = self.QMatoms[i_atom]
              atom.set_xyz( positions_lists[0][i_atom][0] * lengthconv, positions_lists[0][i_atom][1] * lengthconv, positions_lists[0][i_atom][2] * lengthconv )

          # now set positions in self.QMother ..
          for i_atom in range(len(self.QMother)):
              atom = self.QMother[i_atom]
              atom.set_xyz( positions_lists[1][i_atom][0] * lengthconv, positions_lists[1][i_atom][1] * lengthconv, positions_lists[1][i_atom][2] * lengthconv )

          # finally set positions in self.QMdrude ..
          for i_atom in range(len(self.QMdrude)):
              atom = self.QMdrude[i_atom]
              atom.set_xyz( positions_lists[2][i_atom][0] * lengthconv, positions_lists[2][i_atom][1] * lengthconv, positions_lists[2][i_atom][2] * lengthconv )


    #*************************************
    # here we setup the Psi4 geometry class
    # geomlist input is a 2D list ( natoms x 4 ) with type , x , y , z for every atom
    def set_geometry( self ):
          print("Setting charge and spin in QM calculations : " , self.QMcharge , self.QMspin )
          if self.QMspin > 1:
             # set UKS
             core.set_local_option('SCF','REFERENCE','UKS')

          #*************** Add MM charges in QMregion for analytic embedding
          # these atoms are in QMother, which contains QMregion minus QMatoms
          Chrgfield = QMMM()  # this class should be available with 'from psi4.driver import *'
          
          for atom in self.QMother:
              Chrgfield.extern.addCharge( atom.charge , atom.x , atom.y , atom.z )   
          core.set_global_option_python('EXTERN', Chrgfield.extern)
 
          # construct geometry string
          geometrystring=' \n '
          geometrystring= geometrystring + str(self.QMcharge) + " " + str(self.QMspin) +" \n"
          # don't reorient molecule.  Don't want symmetry anyway for DMA
          geometrystring= geometrystring + " noreorient  \n  " + " nocom  \n  "

          #**************** Input QMatoms/positions for QM calculation
          for atom in self.QMatoms:
              geometrystring = geometrystring + " " + str(atom.element) + " " + str(atom.x) + " " + str(atom.y) + " " + str(atom.z) + " \n"
          geometrystring = geometrystring + ' symmetry c1 \n '
          # now create Psi4 geometry object
          self.geometry = psi4.geometry( geometrystring )


    # this function computes the DMA from the QM wavefunction, and uses mpfit to get charges on QM atoms
    def get_DMA_charges( self ):
          # print formatted checkpoint file with wavefunction info 
          fw = core.FCHKWriter(self.wfn)
          fchkfile = 'molecule' + '.fchk'
          fw.write(fchkfile)
          p = subprocess.Popen("./fit_charges.sh", stdout=subprocess.PIPE)
          out, err = p.communicate()
          charges=out.splitlines()
          # make sure charges were fit, wavefunction converged, etc.
          # this isn't a perfect check, but it makes sure we have something...
          if len(charges) == len(self.atoms_type): 
                self.charges=[ float(charges[j].decode("utf-8")) for j in range(len(charges)) ]

    #*******************************************
    # this is a wrapper that calls numpy routines to
    # create PME grid positions in Bohr.
    # input is the number of grid points
    # and the box vectors (a list of lists)
    #*******************************************
    def set_PMEgrid_xyz( self , box ):
          # getting increments in every dimension
          xs = np.linspace( 0 , sum([box[0][i]**2 for i in range(3)])**0.5 , self.pme_grid_size , endpoint=False )
          ys = np.linspace( 0 , sum([box[1][i]**2 for i in range(3)])**0.5 , self.pme_grid_size , endpoint=False )
          zs = np.linspace( 0 , sum([box[2][i]**2 for i in range(3)])**0.5 , self.pme_grid_size , endpoint=False )
          # collecting meshgrid arrays
          X, Y, Z = np.meshgrid( xs , ys , zs , indexing='ij' )
          # flattening meshgrid arrays
          x = X.flatten()[:,np.newaxis]
          y = Y.flatten()[:,np.newaxis]
          z = Z.flatten()[:,np.newaxis]

          #pos_lists = []
          #for i in range(pme_grid_size):
          #     for j in range(pme_grid_size):
          #          for k in range(pme_grid_size):
          #               index = [i,j,k]
          #               pos_list = [0,0,0]
          #               for n in range(3):
          #                    for m in range(3):
          #                         pos_list[m] += ( index[n] / pme_grid_size ) * box[n][m]
          #               pos_lists.append(pos_list)
          #self.pmegrid_xyz = np.array( pos_lists )

          #setting results
          self.pmegrid_xyz = np.concatenate( (x,y,z) , axis=1 )
          self.box = box

    #*******************************************
    # this method collects the PME grid points to exclude
    # from the vext claculation in reciprocal space.
    # input is the DFT functional method (a string)
    #*******************************************
    def get_PMEexclude( self ):
          # collecting the DFT quadrature grid
          sup_func = dft.build_superfunctional( self.DFT_functional , True )[0]
          basis = core.BasisSet.build( self.geometry , "ORBITAL" , self.basis )
          Vpot = core.VBase.build( basis , sup_func , "RV" )
          Vpot.initialize()
          x, y, z, w = Vpot.get_np_xyzw()
          quadrature_grid=[]
          for i in range(len(x)):
                quadrature_grid.append( [ x[i] , y[i] , z[i] ] )
          quadrature_grid = np.array( quadrature_grid )
          self.quadrature_grid = quadrature_grid
          # projecting box to reciprocal space
          inverse_box = calculate_inverse_box_vectors( self.box )
          # projecting quadrature grid to reciprocal space
          quadrature_grid_project = project_to_PME_grid( quadrature_grid , inverse_box , self.pme_grid_size )
          xi = quadrature_grid_project
          # performing floor and ceil operations for each dimension
          xx = [np.floor(xi[:,0].reshape((-1,1))),np.ceil(xi[:,0].reshape((-1,1)))]
          xy = [np.floor(xi[:,1].reshape((-1,1))),np.ceil(xi[:,1].reshape((-1,1)))]
          xz = [np.floor(xi[:,2].reshape((-1,1))),np.ceil(xi[:,2].reshape((-1,1)))]
          # populating with every permutation of the above operations
          xs = []
          for i in range(2):
                for j in range(2):
                      for k in range(2):
                            xs.append(np.concatenate((xx[i],xy[j],xz[k]),axis=1))
          # removing non-unique points
          xf = np.unique(np.concatenate(tuple(xs),axis=0),axis=0)
          #xf = np.random.randint(1,150,size=(700,3))
          xf[xf==self.pme_grid_size] = 0
          self.PMEexclude_list = xf

    #*******************************************
    # this method alters the vext grid in order to 
    # account for exclusions.
    # input is the pme_alpha value in nm^-1, nm to 
    # bohr conversion factor, the real positions 
    # of the atoms, and the original vext grid
    #*******************************************
    def set_PMEexclude( self , real_pos , vext ):
          # conversion factor for nm to Bohr
          vext = np.array( vext ) / hartree_to_kjmol
          # preparing an array of indices which will correspond to Vext_correction and casting to type int
          indices = (self.PMEexclude_list[:,0]*self.pme_grid_size**2 + self.PMEexclude_list[:,1]*self.pme_grid_size + self.PMEexclude_list[:,2]).astype(np.int)
          # collecting the positions and charges of atoms in the QMregion
          positions = []
          charges = []

          for i , atom in enumerate(self.QMatoms):
                positions.append(real_pos[0][i])
                charges.append(atom.charge)
          for i , atom in enumerate(self.QMother):
                positions.append(real_pos[1][i])
                charges.append(atom.charge)
          for i , atom in enumerate(self.QMdrude):
                positions.append(real_pos[2][i])
                charges.append(atom.charge)

          # casting positions to n_atoms x 1 x 3 array, which will allow the pmegrid positions in realspace to be broadcast onto the QMregion positions later
          positions = np.array(positions)[:,np.newaxis,:] * nm_to_bohr
          charges = np.array(charges)
          # getting realspace pmegrid positions that are in QMregion
          pme_grid = self.pmegrid_xyz[indices,:]
          # getting least mirror postions between pmegrid positions and the QMregion atom positions, which will broadcast to produce an n_atom x n_gridpoints x 3 array
          dr = get_least_mirror_pos( pme_grid , positions , self.box )
          #dr = positions - pme_grid
          # getting inverse distance, which will produce an n_atom x n_gridpoint array
          inv_dr = 1 / np.linalg.norm( dr , axis=2 )
          pme_alpha_bohr = self.pme_alpha / ( nm_to_bohr )
          alpha_dr = pme_alpha_bohr / ( inv_dr )
          # subtracting contributions from the indexed Vext_correction as an n_gridpoint x 1 array
          vext[indices] -= np.sum( (charges[np.newaxis,:]*inv_dr.T*erf(alpha_dr).T) , axis=1 ) 

          # adding correction for switching from gaussians to point charges at the boundary of the QM and MM regions
          if self.pme_real_correction:
                vext[indices] += np.sum( (charges[np.newaxis,:]*inv_dr.T*erfc(alpha_dr).T) , axis=1 )

          return vext

    #*******************************************
    # this method performs the Psi4 Energy
    # energy calculation call
    #*******************************************
    def calc_energy( self , vext=None , box=None ):

          if self.qmmm_ewald:
              core.set_local_option("SCF","QMMM", self.qmmm_ewald)
              ( self.energy , self.wfn ) = psi4.energy( self.DFT_functional , return_wfn=True , pme_grid_size=self.pme_grid_size , vexternal_grid=vext , box=box , interpolation_method="interpn" )

          else:
              core.set_local_option("SCF","QMMM", self.qmmm_ewald)
              ( self.energy , self.wfn ) = psi4.energy( self.DFT_functional , return_wfn=True )

#****************************************************
# this is a standalone helper method, outside of class.
# input is two coordinate arrays that may be broadcast 
# together and the box vectors (a list of lists)
#****************************************************
def get_least_mirror_pos( i_vec, j_vec, box_vec ):
    # broadcasting distance between elements of the two arrays
    r = i_vec - j_vec
    r -= (np.array(box_vec[2])[:,np.newaxis,np.newaxis]*np.floor(r[:,:,2]/box_vec[2][2]+0.5).T).T
    r -= (np.array(box_vec[1])[:,np.newaxis,np.newaxis]*np.floor(r[:,:,1]/box_vec[1][1]+0.5).T).T
    r -= (np.array(box_vec[0])[:,np.newaxis,np.newaxis]*np.floor(r[:,:,0]/box_vec[0][0]+0.5).T).T

    return r

# this is standalone helper method, outside of class
# input is an OpenMM state object
def get_box_vectors_Bohr( state ):
    # this is returned as a list of Openmm Vec3 filled with quantities (value, unit)
    # convert to unitless list to pass to Psi4
    box_temp = state.getPeriodicBoxVectors()
    box=[]
    for i in range(3):
        box.append( [ box_temp[i][0]._value * nm_to_bohr , box_temp[i][1]._value * nm_to_bohr, box_temp[i][2]._value * nm_to_bohr ] )

    return box

