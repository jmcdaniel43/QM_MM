from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#******** exclusions for force field 
from shared.MM_exclusions_base import *


#*************************************************
# This is the base MM parent class that is meant for general use when invoking OpenMM
#
# Any specialized/simulation-specific run-control should be implemented in a child class
# of this parent class.  Because this parent class may be used for different types of
# simulations, it may be utilized in several different github repositories which may
# be best managed using Git subtrees
#
# We set reasonable default run-parameters in this base class, and allow
# modification to these arguments with **kwargs input.
# 
#**************************************************
class MM_base(object):
    # required input: 1) list of pdb files, 2) list of residue xml files, 3) list of force field xml files.
    def __init__( self , pdb_list , residue_xml_list , ff_xml_list , **kwargs  ):
        #*************************************
        #  DEFAULT RUN PARAMETERS: input in **kwargs may overide defaults
        #**************************************
        self.temperature = 300*kelvin
        self.temperature_drude = 1*kelvin
        self.friction = 1/picosecond
        self.friction_drude = 1/picosecond
        self.timestep = 0.001*picoseconds
        self.small_threshold = 1e-6  # threshold for charge magnitude
        self.cutoff = 1.4*nanometer

        # reading inputs from **kwargs
        if 'temperature' in kwargs :
            self.temperature = int(kwargs['temperature'])*kelvin
        if 'temperature_drude' in kwargs :
            self.temperature_drude = int(kwargs['temperature_drude'])*kelvin
        if 'friction' in kwargs :
            self.friction = int(kwargs['friction'])/picosecond
        if 'friction_drude' in kwargs :
            self.friction_drude = int(kwargs['friction_drude'])/picosecond
        if 'timestep' in kwargs :
            self.timestep = float(kwargs['timestep'])*picoseconds
        if 'small_threshold' in kwargs :
            self.small_threshold = float(kwargs['small_threshold'])
        if 'cutoff' in kwargs :
            self.cutoff = float(kwargs['cutoff'])*nanometer


        # load bond definitions before creating pdb object (which calls createStandardBonds() internally upon __init__).  Note that loadBondDefinitions is a static method
        # of Topology, so even though PDBFile creates its own topology object, these bond definitions will be applied...
        for residue_file in residue_xml_list:
            Topology().loadBondDefinitions(residue_file)

        # now create pdb object, use first pdb file input
        self.pdb = PDBFile( pdb_list[0] )

        # create modeller
        self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
        # create force field
        self.forcefield = ForceField(*ff_xml_list)
        # add extra particles
        self.modeller.addExtraParticles(self.forcefield)

        # If QM/MM, add QMregion to topology for exclusion in vext calculation...
        #if self.QMMM :
        #    self.modeller.topology.addQMatoms( self.QMregion_list )

        # polarizable simulation?  Figure this out by seeing if we've added any Drude particles ...
        self.polarization = True
        if self.pdb.topology.getNumAtoms() == self.modeller.topology.getNumAtoms():
            self.polarization = False

        if self.polarization :
            #************** Polarizable simulation, use Drude integrator with standard settings
            self.integrator = DrudeLangevinIntegrator(self.temperature, self.friction, self.temperature_drude, self.friction_drude, self.timestep)
            # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
            self.integrator.setMaxDrudeDistance(0.02)
        else :
            #************** Non-polarizable simulation
            self.integrator = LangevinIntegrator(self.temperature, self.friction, self.timestep)


        # create openMM system object
        self.system = self.forcefield.createSystem(self.modeller.topology, nonbondedCutoff=self.cutoff, constraints=HBonds, rigidWater=True)
        # get force types and set method
        self.nbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == NonbondedForce][0]
        self.customNonbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomNonbondedForce][0]
        if self.polarization :
            self.drudeForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == DrudeForce][0]
            # will only have this for certain molecules
            self.custombond = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomBondForce][0]

        # set long-range interaction method
        self.nbondedForce.setNonbondedMethod(NonbondedForce.PME)
        self.customNonbondedForce.setNonbondedMethod(min(self.nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))



    #*********************************************
    # set output frequency for coordinate dcd file
    #*********************************************
    def set_trajectory_output( self, filename , write_frequency ):
        self.simmd.reporters = []
        self.simmd.reporters.append(DCDReporter(filename, write_frequency))



    #*********************************************
    # this sets the force groups to be used with PBC
    # call this if a molecule/residue is broken up over PBC,
    # e.g. graphene electrode ...
    #*********************************************
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


    #*********************************************
    # set the platform/OpenMM kernel and initialize simulation object
    #*********************************************
    #*********** Currently can only use 'Reference' for QM/MM ...
    def set_platform( self, platformname ):
        if platformname == 'Reference':
            self.platform = Platform.getPlatformByName('Reference')
            self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
        elif platformname == 'CPU':
            self.platform = Platform.getPlatformByName('CPU')
            self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
        elif platformname == 'OpenCL':
            self.platform = Platform.getPlatformByName('OpenCL')
            # we found weird bug with 'mixed' precision on OpenCL related to updating parameters in context for gold/water simulation...
            #self.properties = {'OpenCLPrecision': 'mixed'}
            self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform)
        elif platformname == 'CUDA':
            self.platform = Platform.getPlatformByName('CUDA')
            self.properties = {'Precision': 'mixed'}
            self.simmd = Simulation(self.modeller.topology, self.system, self.integrator, self.platform, self.properties)
        else:
            print(' Could not recognize platform selection ... ')
            sys.exit(0)
        self.simmd.context.setPositions(self.modeller.positions)



    #***************************************
    # this generates force field exclusions that we commonly utilize for water simulations
    #
    # if flag_SAPT_FF_exclusions=True, then will also set exclusions for SAPT-FF force field...
    #***************************************
    def generate_exclusions(self, water_name = 'HOH', flag_hybrid_water_model = False ,  flag_SAPT_FF_exclusions = True ):

        # if special exclusion for SAPT-FF force field ...
        if flag_SAPT_FF_exclusions:
            generate_SAPT_FF_exclusions( self )

        # if using a hybrid water model, need to create interaction groups for customnonbonded force....
        if flag_hybrid_water_model:
            generate_exclusions_water(self.simmd, self.customNonbondedForce, water_name )

        # having both is redundant, as SAPT-FF already creates interaction groups for water/other
        if flag_SAPT_FF_exclusions and flag_hybrid_water_model:
            print( "redundant setting of flag_SAPT_FF_exclusions and flag_hybrid_water_model")
            sys.exit()


        # now reinitialize to make sure changes are stored in context
        state = self.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
        positions = state.getPositions()
        self.simmd.context.reinitialize()
        self.simmd.context.setPositions(positions)



