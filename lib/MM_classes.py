from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import *
#******** contains parent class
from shared.MM_class_base import *

#*************************************************
# This is a child MM system class for QM/MM simulations
# that inherits from parent MM_base class in MM_class_base.py
#
# IMPORTANT:  must use compiled version of customized OpenMM/Psi4 for certain features ...
#**************************************************
class MM_QMMM(MM_base):
    # required input: 1) list of pdb files, 2) list of residue xml files, 3) list of force field xml files.
    def __init__( self , pdb_list , residue_xml_list , ff_xml_list , **kwargs  ):
        #storing inputs for later
        self.inputs = kwargs
        self.qmmm_ewald = False
        
        # constructor for Parent...
        super().__init__( pdb_list , residue_xml_list , ff_xml_list , **kwargs )


        # inputs from **kwargs
        if 'QMatoms_range' in kwargs :
            ( QMatoms_low , QMatoms_high ) = kwargs['QMatoms_range'].split(',')
            self.QMatoms_list = range( int( QMatoms_low.strip( ' (),' ) ) , int( QMatoms_high.strip( ' (),' ) ) )
        if 'qmmm_ewald' in kwargs :
            self.qmmm_ewald = eval(kwargs['qmmm_ewald'])
        if 'qmmm_cutoff' in kwargs :
            self.qmmm_cutoff = float(kwargs['qmmm_cutoff'])

        if self.qmmm_ewald :
            if 'pme_grid_size' in kwargs :
                self.pme_grid_size = int(kwargs['pme_grid_size'])
            if 'pme_alpha' in kwargs :
                self.pme_alpha = float(kwargs['pme_alpha'])


        # this sets the PME parameters in OpenMM.  The grid size is important for the accuracy of the 
        # external potental in the DFT quadrature, since this is interpolated from the PME grid
        if self.qmmm_ewald:
            self.nbondedForce.setPMEParameters( self.pme_alpha , self.pme_grid_size , self.pme_grid_size , self.pme_grid_size )



    #*********************************************
    # set the platform/OpenMM kernel and initialize simulation object
    # this overides method in parent class
    # with options for QM/MM/PME in properties
    #*********************************************
    def set_platform( self, platformname ):
        if platformname == 'Reference':
            self.platform = Platform.getPlatformByName('Reference')
            if self.qmmm_ewald :
                self.properties = {'ReferenceVextGrid': 'true'}
                self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
            else :
                self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
        elif platformname == 'CPU':
            self.platform = Platform.getPlatformByName('CPU')
            if self.qmmm_ewald :
                self.properties = {'ReferenceVextGrid': 'true'}
                self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
            else :
                self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
        elif platformname == 'OpenCL':
            self.platform = Platform.getPlatformByName('OpenCL')
            if self.qmmm_ewald :
                print( 'Can only run QM/MM simulation with reference/CPU platforms !')
                sys.exit()
            else :
                # we found weird bug with 'mixed' precision on OpenCL related to updating parameters in context for gold/water simulation...
                #self.properties = {'OpenCLPrecision': 'mixed'}
                self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
        elif platformname == 'CUDA':
            self.platform = Platform.getPlatformByName('CUDA')
            self.properties = {'Precision': 'mixed'}
            if self.qmmm_ewald :
                print( 'Can only run QM/MM simulation with reference/CPU platforms !')
                sys.exit()
            else :
                self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
        else:
            print(' Could not recognize platform selection ... ')
            sys.exit(0)
        self.simmd.context.setPositions(self.modeller.positions)



    #**********************************
    # Reinitializes the QM system object
    # kwargs is a dictionary of arguments
    #**********************************
    def reset( self , **kwargs ):
          self.__init__( **kwargs )

    def set_QMregion_parameters( self , QMatoms_list , qmmm_cutoff ):
          self.qmmm_cutoff = qmmm_cutoff
          self.QMatoms_list = QMatoms_list



    #**************************
    # input atom_lists is list of lists of atoms
    # returns a list of lists of elements, charges with one-to-one correspondence...
    #**************************
    def get_element_charge_for_atom_lists( self, atom_lists ):

        element_lists=[]
        charge_lists=[]
        # loop over lists in atom_lists , and add list to element_lists , charge_lists
        for atom_list in atom_lists:
            element_list=[]
            charge_list=[]
            # loop over atoms in topology and match atoms from list...
            for atom in self.simmd.topology.atoms():
                # if in atom_list ..
                if atom.index in atom_list:
                    element = atom.element
                    # get atomic charge from force field...
                    (q_i, sig, eps) = self.nbondedForce.getParticleParameters(atom.index)
                    # add to lists
                    if len(dir(element)) == 40:
                         element_list.append( element.symbol )
                    else:
                         element_list.append( atom.name )
                    charge_list.append( q_i._value )

            # now add to element_lists , charge_lists ..
            element_lists.append( element_list )
            charge_lists.append( charge_list )

        return element_lists , charge_lists

  
    #**************************
    # input atom_lists is list of lists of atoms
    # returns a list of lists of positions with one-to-one correspondence...
    #**************************
    def get_positions_for_atom_lists( self , atom_lists ):

        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
        positions = state.getPositions()
        QM_pos = positions[self.QMatoms_list[0]]._value
        box_vectors = [state.getPeriodicBoxVectors()[i]._value for i in range(3)]
        lm_position_lists=[]
        real_position_lists=[]
        # loop over lists in atom_lists , and add list to position_lists
        for atom_list in atom_lists:
            lm_position_list=[]
            real_position_list=[]
            for index in atom_list:
                ind_pos = positions[index]._value
                real_position_list.append( ind_pos )
                r = get_least_mirror_pos( ind_pos , QM_pos , box_vectors)
                lm_position_list.append( [r[i]+QM_pos[i] for i in range(3)] )
            # now add to position_lists ...
            lm_position_lists.append( lm_position_list )
            real_position_lists.append( real_position_list )

        return lm_position_lists , real_position_lists



    #**************************
    # this method gets the list of atoms that are in the QM region, including the specified QM system atoms
    # input is a tuple of atom indices for the QM system and a cutoff distance in nanometers
    #**************************
    def get_QMregion_list(self):
          # converting Angstrom cutoff to nm
          cutoff = self.qmmm_cutoff/10.0
          # getting the atom index for the first atom listed for each residue
          res_atom_ind = []
          res_list = [res for res in self.simmd.topology.residues()]
          for res in res_list:
               res_atom_ind.append(res._atoms[0].index)
          # getting current box size for minimum mirror image calculation
          state = self.simmd.context.getState( getEnergy=False , getForces=False , getVelocities=True , getPositions=True , getParameters=True )
          pos = state.getPositions()
          box_vectors = [state.getPeriodicBoxVectors()[j]._value for j in range(3)]
          QM_pos = [sum( [pos[i]._value[j] for i in self.QMatoms_list] ) / len( self.QMatoms_list ) for j in range(3)]
          # populating QMregion_list and QMdrudes_list
          QMregion_list = []
          QMdrudes_list = []
          for i in range(len(pos)):
               r = get_least_mirror_pos(QM_pos,pos[i]._value,box_vectors)
               dist = sum([r[i]**2 for i in range(3)])**(0.5)
               if dist < cutoff and i in res_atom_ind:
                     ind = res_atom_ind.index(i)
                     if i in self.QMatoms_list:
                         QM_atoms = [atom.index for atom in res_list[ind]._atoms]
                         QM_elements = [atom.element for atom in res_list[ind]._atoms]
                         for j in range(len(QM_atoms)):
                             if hasattr(QM_elements[j],'symbol'):
                                 QMregion_list.append(QM_atoms[j])
                             else:
                                 QMdrudes_list.append(QM_atoms[j])
                     else:
                          QMregion_list.extend([atom.index for atom in res_list[ind]._atoms])

          return tuple( QMregion_list ) , tuple( QMdrudes_list )


    #********************
    # this method writes a pdb file
    # using the current state
    #*******************
    def write_pdb( self , name ):
          state = self.simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,enforcePeriodicBox=True)
          pos_list = state.getPositions()
          self.simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
          PDBFile.writeFile( self.simmd.topology , pos_list , open( self.out_dir.split( '.' )[0] + '_' + name + '.pdb', 'w' ) )

#****************************************************
# this is a standalone helper method, outside of class
# input is three list vectors
#****************************************************
def get_least_mirror_pos( i_vec, j_vec, box_vec ):
    # getting index of QM atoms in the QMregion_list to avoid self interaction calculation
    r = i_vec - j_vec
    r -= box_vec[2]*math.floor(r[2]/box_vec[2][2]+0.5)
    r -= box_vec[1]*math.floor(r[1]/box_vec[1][1]+0.5)
    r -= box_vec[0]*math.floor(r[0]/box_vec[0][0]+0.5)

    return r
