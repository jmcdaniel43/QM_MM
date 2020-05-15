from sys import stdout
import psi4
import psi4.core as core
from psi4.driver import *

import subprocess

# This Module defines new classes that we use to run a QM/MM simulation

#*************************** README  **************************************
#  This module defines QM classe that is used in a QM/MM simulation
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
    def __init__(self, QMname, basis, dft_spherical_points , dft_radial_points , scf_type , qmmm='true' ):
          # setup options for DMA
          radii=['H',0.53,'C',0.53,'N',0.53,'O',0.53]
          core.set_local_option('GDMA', 'GDMA_RADIUS', radii )
          core.set_local_option('GDMA', 'GDMA_LIMIT', 4 )
          core.set_local_option('GDMA', 'GDMA_MULTIPOLE_UNITS', 'AU' )
          core.set_local_option('GDMA', 'GDMA_SWITCH', 0.0 )
          # name
          self.molname = QMname

          # set basis , quadrature settings, scf_type, qmmm ...
          psi4.set_options( {'basis' : basis , 'dft_spherical_points' : dft_spherical_points , 'dft_radial_points' : dft_radial_points , 'scf_type' : scf_type , 'qmmm' : qmmm  } )


    #**********************************
    # QMatoms_list is tuple specifying atom indices of QM atom
    # QMother_list specifies additional MM embedding atoms
    # element_lists / charge_lists correspond to [ QMatoms_list , QMother_list ] ...
    #**********************************
    def set_QM_region(self, element_lists , charge_lists , QMatoms_list , QMother_list ):
          # create lists of atom_QM objects: there are two lists here:
          # QMatoms stores data for QMatoms
          # QMother stores data for MM embedding atoms
          self.QMatoms = []
          self.QMother = []

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



    #*************************************
    # here we setup the Psi4 geometry class
    # geomlist input is a 2D list ( natoms x 4 ) with type , x , y , z for every atom
    def set_geometry(self, charge, spin ):
          print("Setting charge and spin in QM calculations : " , charge , spin )
          if spin > 1:
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
          geometrystring= geometrystring + str(charge) + " " + str(spin) +" \n"
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




# this is standalone helper method, outside of class
# input is an OpenMM state object
def get_box_vectors_Bohr( state, nm_to_bohr ):
    # this is returned as a list of Openmm Vec3 filled with quantities (value, unit)
    # convert to unitless list to pass to Psi4
    box_temp = state.getPeriodicBoxVectors()
    box=[]
    for i in range(3):
        box.append( [ box_temp[i][0]._value * nm_to_bohr , box_temp[i][1]._value * nm_to_bohr, box_temp[i][2]._value * nm_to_bohr ] )

    return box

