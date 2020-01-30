from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import psi4
import psi4.core as core
from psi4.driver import *

import subprocess

# This Module defines new classes that we use to run a QM/MM simulation

#*************************** README  **************************************
#  This module defines classes that are used in a QM/MM simulation
#  Interfacing Psi4  and OpenMM
#
#  Because these codes use different units and datastructures, conversions
#  between these are done within the QM/MM classes, so that input and output
#  passed to/from the classes is given in units/datastructures of the
#  particular libraries (Psi4, OpenMM) that it comes from/goes to. 
#
#**************************************************************************

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


    # Input QMregion is tuple specifying atom indices of QM region, which
    # encompasses both QMatoms (input tuple) and additional MM embedding atoms
    # MMsys is input to figure out element types of QM atoms
    def set_QM_region(self, MMsys , QMregion , QMatoms ):
          # store data structures for later...
          self.QMregion = QMregion
          self.QMatoms = QMatoms
          self.atoms_index=[]
          self.atoms_type=[]
          self.atoms_charge=[]
          # maybe not the best way to do it, but loop over all MM atoms ...
          for atom in MMsys.simmd.topology.atoms():
             # if in QM region ..
             if atom.index in QMregion:
                 element = atom.element
                 # get partial charge from force field.  Need this if this is MM atom in QMregion...
                 (q_i, sig, eps) = MMsys.nbondedForce.getParticleParameters(atom.index)

                 self.atoms_charge.append( q_i._value )
                 self.atoms_index.append( atom.index )
                 self.atoms_type.append( element.symbol )


       # here we setup the Psi4 geometry class
       # geomlist input is a 2D list ( natoms x 4 ) with type , x , y , z for every atom
       # input is in OpenMM units, and internally converted to Psi4 units
    def set_geometry(self, geomlist , charge, spin ):
          lengthconv = 10.0  # nm to angstrom
          print("Setting charge and spin in QM calculations : " , charge , spin )
          if spin > 1:
             # set UKS
             core.set_local_option('SCF','REFERENCE','UKS')

          #*************** Add MM charges in QMregion for analytic embedding
          # MMembed atoms should start after QMatoms in QMregion ...
          Chrgfield = QMMM()  # this class should be available with 'from psi4.driver import *'
          MMembed_start = len(self.QMatoms)
          for i in range(MMembed_start,len(self.QMregion)):
              Chrgfield.extern.addCharge( self.atoms_charge[i] , geomlist[i][1]*lengthconv , geomlist[i][2]*lengthconv , geomlist[i][3]*lengthconv )  # q , x , y , z     
          core.set_global_option_python('EXTERN', Chrgfield.extern)
 
          # construct geometry string
          geometrystring=' \n '
          geometrystring= geometrystring + str(charge) + " " + str(spin) +" \n"
          # don't reorient molecule.  Don't want symmetry anyway for DMA
          geometrystring= geometrystring + " noreorient  \n  " + " nocom  \n  "
          for i in range(len(self.QMatoms)):
              geometrystring = geometrystring + " " + str(geomlist[i][0]) + " " + str(geomlist[i][1]*lengthconv) + " " + str(geomlist[i][2]*lengthconv) + " " + str(geomlist[i][3]*lengthconv) + " \n"
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



# This class controls the MM region of the simulation
class MM(object):
    # input to init is 3 lists , list of pdb files, list of residue xml files, list of force field xml files , and list of atoms in QMregion   
    def __init__(self, pdb_list , residue_xml_list , ff_xml_list , QMregion , cutoff ):
          # ********************* Standard Simulation settings.  Change these if necessary
          # this controls openMM API kernel
          self.NPT = False
          self.temperature = 300*kelvin
          self.temperature_drude = 1*kelvin
          self.friction = 1/picosecond
          self.friction_drude = 1/picosecond
          self.timestep = 0.001*picoseconds
          self.pressure = Vec3(1.0,1.0,1.0)*atmosphere
          self.cutoff = cutoff 
          self.barofreq = 100
          # set Open MM Integrator, use Drude integrator with standard settings
          self.integrator = DrudeLangevinIntegrator(self.temperature, self.friction, self.temperature_drude, self.friction_drude, self.timestep)    
          # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
          self.integrator.setMaxDrudeDistance(0.02)  
           
          # store QMregion
          self.QMregion = QMregion

          # now create pdb object, use first pdb file input
          self.pdb = PDBFile( pdb_list[0] )
          # update pdb topology with residues
          for residue_file in residue_xml_list:
               self.pdb.topology.loadBondDefinitions(residue_file)
          self.pdb.topology.createStandardBonds()
          # create modeller
          self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
          # create force field
          self.forcefield = ForceField(*ff_xml_list)
          # add extra particles
          self.modeller.addExtraParticles(self.forcefield)

          # Add tuple of QM atoms to topology
          self.modeller.topology.addQMatoms( self.QMregion )

          # create openMM system object
          self.system = self.forcefield.createSystem(self.modeller.topology, nonbondedCutoff=self.cutoff, constraints=None, rigidWater=True)
          # get force types and set method
          self.nbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == NonbondedForce][0]
          self.customNonbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomNonbondedForce][0]
          self.drudeForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == DrudeForce][0]
          self.nbondedForce.setNonbondedMethod(NonbondedForce.PME)
          self.customNonbondedForce.setNonbondedMethod(min(self.nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))
          # add barostat if NPT
          if self.NPT:
               # allow only the z-dimension to change, graphene x/y layer is fixed
               self.barostat = MonteCarloAnisotropicBarostat(self.pressure,self.temperature,False,False,True,self.barofreq)
               self.system.addForce(self.barostat)

    # this sets the force groups to used PBC
    def set_periodic_residue(self, flag):
          for i in range(self.system.getNumForces()):
               f = self.system.getForce(i)
               f.setForceGroup(i)
               # if using PBC
               if flag:
                      # Here we are adding periodic boundaries to intra-molecular interactions.  Note that DrudeForce does not have this attribute, and
                      # so if we want to use thole screening for graphite sheets we might have to implement periodic boundaries for this force type
                      if type(f) == HarmonicBondForce or type(f) == HarmonicAngleForce or type(f) == PeriodicTorsionForce or type(f) == RBTorsionForce:
                            f.setUsesPeriodicBoundaryConditions(True)
                            f.usesPeriodicBoundaryConditions()

    # this sets the PME parameters in OpenMM.  The grid size is important for the accuracy of the external potential
    # in the DFT quadrature, since this is interpolated from the PME grid
    def setPMEParameters( self , pme_alpha , pme_grid_a , pme_grid_b , pme_grid_c ):
        self.nbondedForce.setPMEParameters( pme_alpha , pme_grid_a , pme_grid_b , pme_grid_c )

    # this sets the platform for OpenMM simulation and initializes simulation object
    #*********** Currently can only use 'Reference' for generating vext_grid !
    def set_platform( self, platformname ):
          if platformname == 'Reference':
               self.properties = {'ReferenceVextGrid': 'true'}
               self.platform = Platform.getPlatformByName('Reference')
               self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
          #elif platformname == 'OpenCL':
          #     self.platform = Platform.getPlatformByName('OpenCL')
          #     self.properties = {'OpenCLPrecision': 'mixed','ReferenceVextGrid': 'true'}
          #     self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
          #elif platformname == 'CUDA':
          #     self.platform = Platform.getPlatformByName('CUDA')
          #     self.properties = {'CUDAPrecision': 'mixed','ReferenceVextGrid': 'true'}
          #     self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
          else:
               print(' Only "Reference" platform is currently implemented for QM/MM ... ')
               sys.exit(0)
          self.simmd.context.setPositions(self.modeller.positions)

    # this method sets an umbrella potential constraining the centroid of input molecule "mol1" to distance from atom "atom" on mol2
    # input simmd is an openMM simulation object, input system is an openMM system object, input modeller is openMM modeller object
    def setumbrella(self, mol1, mol2, atomtype, k , r0centroid ):
          # create groups for centroid     
          g1=[]
          g2=[]
          for res in self.simmd.topology.residues():
               if res.name == mol1:
                     for i in range(len(res._atoms)):
                           g1.append(res._atoms[i].index)
                     break
          for res in self.simmd.topology.residues():
               if res.name == mol2:
                     for i in range(len(res._atoms)):
                           if res._atoms[i].name == atomtype:
                                 g2.append(res._atoms[i].index)
                     break
          self.Centroidforce = CustomCentroidBondForce(2,"0.5*k*(distance(g1,g2)-r0centroid)^2") 
          self.system.addForce(self.Centroidforce)
          self.Centroidforce.addPerBondParameter("k")
          self.Centroidforce.addPerBondParameter("r0centroid")
          self.Centroidforce.addGroup(g1)
          self.Centroidforce.addGroup(g2)
          bondgroups =[0,1]
          bondparam = [k,r0centroid]
          self.Centroidforce.addBond(bondgroups,bondparam)
          self.Centroidforce.setUsesPeriodicBoundaryConditions(True)
          self.Centroidforce.addGlobalParameter('r0centroid',r0centroid)
          #self.Centroidforce.addEnergyParameterDerivative('r0centroid')

          for i in range(self.system.getNumForces()):
               f = self.system.getForce(i)
               f.setForceGroup(i)

          self.simmd.context.reinitialize()
          self.simmd.context.setPositions(self.modeller.positions)
          self.simmd.reporters = []
          self.simmd.reporters.append(DCDReporter('md_npt1.dcd', 10))

    # this method gets positions of the QM atoms  
    def get_positions_QM( self, QMsys ):
          # OpenMM state object stores positions
          state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
          # these are all the positions, we only want some...
          position_all = state.getPositions()

          # output a list with QM elements and positions
          position_QM =[]
          for i in range(len(QMsys.atoms_index)):
               index = QMsys.atoms_index[i]
               position_QM.append( [ QMsys.atoms_type[i] , position_all[index][0]._value , position_all[index][1]._value , position_all[index][2]._value ] )
          
          return position_QM   

      
    # this method updates charges on MM atoms
    def update_charges( self , QMsys ):
          # loop over QM atoms, and fill in new charges
          for i in range(len(QMsys.atoms_index)):
               q = QMsys.charges[i]
               index = QMsys.atoms_index[i]
               # for SAPT-FF, set sigma=1, epsilon=0 (no LJ interaction)
               self.nbondedForce.setParticleParameters(index, q, 1.0 , 0.0)
          # now update the parameters in the context object
          self.nbondedForce.updateParametersInContext( self.simmd.context )


