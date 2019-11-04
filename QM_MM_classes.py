from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import psi4
import psi4.core as core
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
    def __init__(self, name):
          # setup options for DMA
          radii=['H',0.53,'C',0.53,'N',0.53,'O',0.53]
          core.set_local_option('GDMA', 'GDMA_RADIUS', radii )
          core.set_local_option('GDMA', 'GDMA_LIMIT', 4 )
          core.set_local_option('GDMA', 'GDMA_MULTIPOLE_UNITS', 'AU' )
          core.set_local_option('GDMA', 'GDMA_SWITCH', 0.0 )
          # initialize External Field
          core.set_local_option('SCF','PERTURB_H',1)
          core.set_local_option('SCF','PERTURB_WITH','DIPOLE')
          # name
          self.molname = name

    # this creates atomlist that defines QM region.  Uses MM object to find atoms
    def set_QM_region(self, MMsys ):
          self.atoms_index=[]
          self.atoms_type=[]
          for res in MMsys.simmd.topology.residues():
             if res.name == self.molname:
                 for i in range(len(res._atoms)):
                 # if not a Drude oscillator, add to QM
                     if 'D' in res._atoms[i].name:
                         continue
                     else:
                         self.atoms_index.append( res._atoms[i].index )
                         element = res._atoms[i].element
                         self.atoms_type.append( element.symbol )
                 break


       # here we setup the Psi4 geometry class
       # geomlist input is a 2D list ( natoms x 4 ) with type , x , y , z for every atom
       # input is in OpenMM units, and internally converted to Psi4 units
    def set_geometry(self, geomlist , charge ):
          lengthconv = 10.0  # nm to angstrom
          # figure out spin based on charge, allow 0, -1, 1       
          if charge == 0:
             spin=1
          elif charge == 1:
             spin=2
             # set UKS
             core.set_local_option('SCF','REFERENCE','UKS')
          elif charge == -1:
             spin=2
             # set UKS
             core.set_local_option('SCF','REFERENCE','UKS')
          else:
             print( 'unrecognized charge setting in QM.geometry()' )
             sys.exit(1)
          # construct geometry string
          geometrystring=' \n '
          geometrystring= geometrystring + str(charge) + " " + str(spin) +" \n"
          # don't reorient molecule.  Don't want symmetry anyway for DMA
          geometrystring= geometrystring + " noreorient  \n  " + " nocom  \n  "
          for i in range(len(geomlist)):
              geometrystring = geometrystring + " " + str(geomlist[i][0]) + " " + str(geomlist[i][1]*lengthconv) + " " + str(geomlist[i][2]*lengthconv) + " " + str(geomlist[i][3]*lengthconv) + " \n"
          geometrystring = geometrystring + ' symmetry c1 \n '
          # now create Psi4 geometry object
          self.geometry = psi4.geometry( geometrystring )

       # set external field for QM
    def set_external_field( self , Ex , Ey, Ez ):
          # kj/mol/e/nm from openMM to au for Psi4
          fieldconv = 1.0 / 2625.5 / 10 / 1.8897 ;
          core.set_local_option('SCF','PERTURB_DIPOLE',[ Ex*fieldconv , Ey*fieldconv , Ez*fieldconv ] )
 


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
    # input to init is 4 lists, list of pdb files, list of residue xml files, list of force field xml files, and list of force field xml files for Efield calculation     
    def __init__(self, pdb_list , residue_xml_list , ff_xml_list, ff_xml_Efield_list):
          # ********************* Standard Simulation settings.  Change these if necessary
          # this controls openMM API kernel
          self.NPT = True
          self.temperature = 300*kelvin
          self.temperature_drude = 1*kelvin
          self.friction = 1/picosecond
          self.friction_drude = 1/picosecond
          self.timestep = 0.001*picoseconds
          self.pressure = Vec3(1.0,1.0,1.0)*atmosphere
          self.cutoff = 1.4*nanometer  
          self.barofreq = 100
          # set Open MM Integrator, use Drude integrator with standard settings
          self.integrator = DrudeLangevinIntegrator(self.temperature, self.friction, self.temperature_drude, self.friction_drude, self.timestep)    
          # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
          self.integrator.setMaxDrudeDistance(0.02)  
          # this integrator won't be used, its just for computing external electric field
          self.integrator_junk = DrudeLangevinIntegrator(self.temperature, self.friction, self.temperature_drude, self.friction_drude, self.timestep)
           
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
          # create force field for Electric Field
          self.forcefield2 = ForceField(*ff_xml_Efield_list)
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

          # now setup system2 for Efield
          self.system2 = self.forcefield2.createSystem(self.modeller.topology, nonbondedCutoff=self.cutoff, constraints=None, rigidWater=True)
          self.nbondedForce2 = [f for f in [self.system2.getForce(i) for i in range(self.system2.getNumForces())] if type(f) == NonbondedForce][0]
          self.nbondedForce2.setNonbondedMethod(NonbondedForce.PME)


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

    # this sets the platform for OpenMM simulation and initializes simulation object
    def set_platform( self, platformname ):
          if platformname == 'CPU':
               self.platform = Platform.getPlatformByName('CPU')
               self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
               self.simfield = Simulation(self.modeller.topology, self.system2, self.integrator_junk, self.platform)
          elif platformname == 'OpenCL':
               self.platform = Platform.getPlatformByName('OpenCL')
               self.properties = {'OpenCLPrecision': 'mixed','OpenCLDeviceIndex':'0'}
               self.properties2 = {'OpenCLPrecision': 'mixed','OpenCLDeviceIndex':'1'}
               self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
               self.simfield = Simulation(self.modeller.topology, self.system2, self.integrator_junk, self.platform, self.properties2)
          elif platformname == 'CUDA':
               self.platform = Platform.getPlatformByName('CUDA')
               self.properties = {'CUDAPrecision': 'mixed','CUDADeviceIndex':'0'}
               self.properties2 = {'CUDAPrecision': 'mixed','CUDADeviceIndex':'1'}
               self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
               self.simfield = Simulation(self.modeller.topology, self.system2, self.integrator_junk, self.platform, self.properties2)
          else:
               print(' platform type not recognized in MM.set_platform( self, platformname ) ')
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

    # this method gets positions and field  
    def get_positions_field( self, QMsys ):
          # OpenMM state object stores positions
          state = self.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
          # these are all the positions, we only want some...
          position_all = state.getPositions()
          # store MM energy
          self.energy = state.getPotentialEnergy()

          # output a list with QM elements and positions
          position_QM =[]
          for i in range(len(QMsys.atoms_index)):
               index = QMsys.atoms_index[i]
               position_QM.append( [ QMsys.atoms_type[i] , position_all[index][0]._value , position_all[index][1]._value , position_all[index][2]._value ] )
          
          # now state object for forces (Efield)
          self.simfield.context.setPositions(position_all)
          # forces with Electrostatics only
          state2 = self.simfield.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
          # these are electrostatic forces, so get field
          forces = state2.getForces()

	  # as of now, we can only use constant electric field on QM.
          # choose field on first atom for now
          index = QMsys.atoms_index[0]
          (q, sig, eps) = self.nbondedForce2.getParticleParameters(index)
          field = forces[index] / q

          return position_QM, [ field[0]._value , field[1]._value , field[2]._value ]    
      
    # this method updates charges on MM atoms
    def update_charges( self , QMsys ):
          # loop over QM atoms, and fill in new charges
          for i in range(len(QMsys.atoms_index)):
               q = QMsys.charges[i]
               index = QMsys.atoms_index[i]
               # for SAPT-FF, set sigma=1, epsilon=0 (no LJ interaction)
               self.nbondedForce.setParticleParameters(index, q, 1.0 , 0.0)
               self.nbondedForce2.setParticleParameters(index, q, 1.0 , 0.0) 
          # now update the parameters in the context object
          self.nbondedForce.updateParametersInContext( self.simmd.context )
          self.nbondedForce2.updateParametersInContext( self.simfield.context )


