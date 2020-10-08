from copy import copy, deepcopy
import numpy as np
import sys
from time import gmtime, strftime, time
from datetime import datetime

# important constants
nm_to_bohr = 18.89726
hartree_to_kjmol = 2625.4996

# this subroutine parses the input file and returns a dictionary of arguments that may be passed to the QM and MM classes
#      file_dir           : Input file directory
def parse_input_file( file_dir ):

    # getting file contents
    fh = open( file_dir )
    lines = fh.readlines()
    fh.close()

    # instantiating input argument dictionary
    input_args = {}

    # defining forbidden values
    forbidden_vals = [ '' , ' ' , '\t' , None ]

    # looping through input lines
    for line in lines:

        # disregarding
        if line[0] != '#' and len(line) > 4:

            if '#' in line:

                line = line.split('#')[0]

            ( key , val ) = line.split( ':' )

            # cleaning keys and values
            key = key.strip(' \t\n')
            val = val.strip(' \t\n')

            # dealing with empty inputs
            if val in forbidden_vals:
                input_args[key] = None

            # otherwise, pair keys and values
            else: 
                input_args[key] = val

    # chekcing inputs and returning
    input_args = check_input_args( input_args )

    return input_args

# this subroutine parses the input file and returns a dictionary of arguments that may be passed to the QM and MM classes
#      input_args         : Input argument dictionary
def check_input_args( input_args ):

    # checking args that have reasonable defaults
    if 'temperature' not in input_args.keys() or input_args['temperature'] is None:
        input_args['temperature'] = '300'
        print(' ''temperature'' not specified, using default of 300 Kelvin')
    if 'temperature_drude' not in input_args.keys() or input_args['temperature_drude'] is None:
        input_args['temperature_drude'] = '1'
        print(' ''temperature_drude'' not specified, using default of 1 Kelvin')
    if 'friction' not in input_args.keys() or input_args['friction'] is None:
        input_args['friction'] = '1'
        print(' ''friction'' not specified, using default of 1 per picosecond')
    if 'friction_drude' not in input_args.keys() or input_args['friction_drude'] is None:
        input_args['friction_drude'] = '1'
        print(' ''friction_drude'' not specified, using default of 1 per picosecond')
    if 'timestep' not in input_args.keys() or input_args['timestep'] is None:
        input_args['timestep'] = '0.001'
        print(' ''timestep'' not specified, using default of 0.001 picoseconds')
    if 'small_threshold' not in input_args.keys() or input_args['small_threshold'] is None:
        input_args['small_threshold'] = '1e-6'
        print(' ''small_threshold'' not specified, using default of 1e-6')
    if 'cutoff' not in input_args.keys() or input_args['cutoff'] is None:
        input_args['cutoff'] = '1.4'
        print(' ''cutoff'' not specified, using default of 1.4 nanometers')
    if 'out_dir' not in input_args.keys() or input_args['out_dir'] is None:
        input_args['out_dir'] = 'FV_QMMM_Test.out'
        print(' ''out_dir'' not specified, using default of ''FV_QMMM_Test.out'' ')
    if 'QMcharge' not in input_args.keys() or input_args['QMcharge'] is None:
        input_args['QMcharge'] = '0'
        print(' ''QMcharge'' not specified, using default of 0')
    if 'QMspin' not in input_args.keys() or input_args['QMspin'] is None:
        input_args['QMspin'] = '1'
        print(' ''QMspin'' not specified, using default of 1')
    if 'platform' not in input_args.keys() or input_args['platform'] is None:
        input_args['platform'] = 'CPU'
        print(' ''platform'' not specified, using default of ''CPU'' ')
    if 'collect_charge_data' not in input_args.keys() or input_args['collect_charge_data'] is None:
        input_args['collect_charge_data'] = 'False'
        print(' ''collect_charge_data'' not specified, using default of ''False'' ')
    if 'qmmm_tare' not in input_args.keys() or input_args['qmmm_tare'] is None:
        input_args['qmmm_tare'] = 'False'
        print(' ''qmmm_tare'' not specified, using default of ''False'' ')
    if 'pme_grid_size' not in input_args.keys() or input_args['pme_grid_size'] is None:
        input_args['pme_grid_size'] = '100'
        print(' ''pme_grid_size'' not specified, using default of 100')
    if 'pme_alpha' not in input_args.keys() or input_args['pme_alpha'] is None:
        input_args['pme_alpha'] = '5.0'
        print(' ''pme_alpha'' not specified, using default of 5.0')
    if 'return_system_init_pdb' not in input_args.keys() or input_args['return_system_init_pdb'] is None:
        input_args['return_system_init_pdb'] = 'True'
        print(' ''return_system_init_pdb'' not specified, using default of ''True'' ')
    if 'return_system_final_pdb' not in input_args.keys() or input_args['return_system_final_pdb'] is None:
        input_args['return_system_final_pdb'] = 'True'
        print(' ''return_system_final_pdb'' not specified, using default of ''True'' ')
    if 'write_frequency' not in input_args.keys() or input_args['write_frequency'] is None:
        input_args['write_frequency'] = '1000'
        print(' ''write_frequency'' not specified, using default of 1000')
    if 'simulation_length' not in input_args.keys() or input_args['simulation_length'] is None:
        input_args['simulation_length'] = '0.1'
        print(' ''simulation_length'' not specified, using default of 0.1')

    # checking args that strictly require inputs
    if 'set_periodic_residue' not in input_args.keys() or input_args['set_periodic_residue'] is None:
        print(' ''set_periodic_residue'' not specified, please provide a boolean-type value for this key in the input file')
        sys.exit( )
    if 'qmmm_ewald' not in input_args.keys() or input_args['qmmm_ewald'] is None:
        print(' ''qmmm_ewald'' not specified, please provide a boolean-type value for this key in the input file')
        sys.exit( )
    if 'qmmm_cutoff' not in input_args.keys() or input_args['qmmm_cutoff'] is None:
        print(' ''qmmm_cutoff'' not specified, please provide a float-type value for this key in the input file')
        sys.exit( )
    if 'functional' not in input_args.keys() or input_args['functional'] is None:
        print(' ''functional'' not specified, please provide a string-type value for this key in the input file')
        sys.exit( )
    if 'basis_set' not in input_args.keys() or input_args['basis_set'] is None:
        print(' ''basis_set'' not specified, please provide a string-type value for this key in the input file')
        sys.exit( )
    if 'quadrature_spherical' not in input_args.keys() or input_args['quadrature_spherical'] is None:
        print(' ''quadrature_spherical'' not specified, please provide an int-type value for this key in the input file')
        sys.exit( )
    if 'quadrature_radial' not in input_args.keys() or input_args['quadrature_radial'] is None:
        print(' ''quadrature_radial'' not specified, please provide an int-type value for this key in the input file')
        sys.exit( )
    if 'pdb_dir' not in input_args.keys() or input_args['pdb_dir'] is None:
        print(' ''pdb_dir'' not specified, please provide a string-type value for this key in the input file')
        sys.exit( )
    if 'res_dir' not in input_args.keys() or input_args['res_dir'] is None:
        print(' ''res_dir'' not specified, please provide a string-type value for this key in the input file')
        sys.exit( )
    if 'ff_dir' not in input_args.keys() or input_args['ff_dir'] is None:
        print(' ''ff_dir'' not specified, please provide a string-type value for this key in the input file')
        sys.exit( )
    if 'QMatoms_range' not in input_args.keys() or input_args['QMatoms_range'] is None:
        print(' ''QMatoms_range'' not specified, please provide a 2-element tuple-type value for this key in the input file')
        sys.exit( )

    return input_args

# this subroutine runs the energy calculation, which requires information from both the QM and MM systems
#      QMsys              : Existing QM object
#      MMsys              : Existing MM object
#      tare               : boolean, if true, the conformational energy will be subtracted from the calculation
def run_qmmm( QMsys , MMsys , QMsys_tare = None ):

    # getting QMregion_list using input values
    MMsys.get_QMregion_list()

    # QMother is the difference between lists ..
    QMother_list = np.setdiff1d( np.array( MMsys.QMregion_list ) , np.array( MMsys.QMatoms_list ) )
    QMdrude_list = MMsys.QMdrudes_list

    # get elements/charges of QM region atoms from MMsys ...
    element_lists , charge_lists = MMsys.get_element_charge_for_atom_lists( [ MMsys.QMatoms_list , QMother_list , QMdrude_list ] )

    embed_charge = sum( charge_lists[1] )

    print('QMregion collected')

    # Fill QM region with atoms.
    QMsys.set_QM_region( element_lists , charge_lists , MMsys.QMatoms_list, QMother_list , QMdrude_list )

    # Get QM positions from MMsystem and set them in QMsys object
    positions_lists , real_pos = MMsys.get_positions_for_atom_lists([ MMsys.QMatoms_list , QMother_list , QMdrude_list] )
    QMsys.set_QM_positions( positions_lists )

    print('QM region set')

    # set geometry of QM region
    QMsys.set_geometry( )

    print('Geometry set')

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

    # performing vext corrections for PME-QM/MM
    if MMsys.qmmm_ewald:

        # collecting state information
        state = MMsys.simmd.context.getState( getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True )
        print( 'done calculating energy, force, and vext' )

        # external potential on PME grid
        vext = state.getVext_grid()

        #** box vectors in units of Bohr
        box = get_box_vectors_Bohr( state )

        # setting PME grid positions
        QMsys.set_PMEgrid_xyz( box )

        # getting PME grid points to exclude from reciprocal space
        QMsys.get_PMEexclude( )

        # getting PME grid corrections and adding to vext_tot
        vext = QMsys.set_PMEexclude( real_pos , vext )

        print('Vext set')

        # calculating electronic energy of the QM system
        QMsys.calc_energy( vext , box )

    else:

        # calculating electronic energy of the QM system
        QMsys.calc_energy( )

    print('Initial E calcd')

    # storing electronic energy value
    E = QMsys.energy

    # checking if a QMsys_tare needs to be initialized
    if QMsys_tare is None and QMsys.qmmm_tare:

        # copying input arguments and QMsystem
        tare_input_args = copy( QMsys.inputs )
        QMsys_tare = copy( QMsys )

        # editing inputs for taring and reinitializing copy with edited inputs
        tare_input_args['qmmm_ewald'] = 'False'
        QMsys_tare.reset( **tare_input_args )

    # subtracting conformational energy
    if QMsys_tare is not None:

        # get elements/charges of QM region atoms from MMsys ...
        element_lists = [ element_lists[0] , [] , [] ]
        charge_lists = [ charge_lists[0] , [] , [] ]

        # Fill QM region with atoms.
        QMsys_tare.set_QM_region( element_lists , charge_lists , MMsys.QMatoms_list, [] , [] )

        # Get QM positions from MMsystem and set them in QMsys object
        positions_lists = [ positions_lists[0] , [] , [] ]
        real_pos = [ real_pos[0] , [] , [] ]
        QMsys_tare.set_QM_positions( positions_lists )

        # set geometry of QM region
        QMsys_tare.set_geometry( )

        # calculating conformational electronic energy of the QM system and storing the value
        QMsys_tare.calc_energy( )
        E_base = QMsys_tare.energy

        # calculating the interaction energy
        E -= E_base

    # checking if embedded charge data should be returned
    if MMsys.collect_charge_data:

        return E, embed_charge

    else:

        return E

# this subroutine facilitates the OpenMM Simulation Loop
#      QMsys              : Existing QM object
#      MMsys              : Existing MM object
#      QMsys_tare         : Existing QM object for taring
def run_simulation( QMsys , MMsys , QMsys_tare = None ):

    # checking if a QMsys_tare needs to be initialized
    if QMsys_tare is None and QMsys.qmmm_tare:

        tare_input_args = copy( QMsys.inputs )
        QMsys_tare = copy( QMsys )
        tare_input_args['qmmm_ewald'] = 'False'
        QMsys_tare.reset( **tare_input_args )

    # checking and creating start drudes pdb
    if MMsys.return_system_init_pdb:
        MMsys.write_pdb( 'start_drudes' )

    # getting state to depict initial energies
    state = MMsys.simmd.context.getState( getEnergy=True , getForces=True , getPositions=True ) 

    # printing initial energies
    print( str(state.getKineticEnergy()) )
    print( str(state.getPotentialEnergy()) )

    for j in range( MMsys.system.getNumForces() ):

        f = MMsys.system.getForce( j )
        print( type(f) , str(MMsys.simmd.context.getState( getEnergy=True , groups=2**j ).getPotentialEnergy()) )

    # this needs to be considered
    #MMsys.simmd.reporters = []
    #MMsys.simmd.reporters.append(DCDReporter('md_output.dcd', write_frequency))

    # opening output file
    fh = open( MMsys.out_dir , 'w' )
    fh.write( str(MMsys.qmmm_cutoff) + 'A' )

    print( 'Starting Simulation...' )

    # checking if an MMsys_CPU needs to be initialized
    if MMsys.platformname is not 'CPU' and MMsys.qmmm_ewald:

        CPU_input_args = copy( MMsys.inputs )
        MMsys_CPU = copy( MMsys )
        CPU_input_args['platform'] = 'CPU'
        MMsys_CPU.reset( **CPU_input_args )

    # otherwise, None object
    else:

        MMsys_CPU = None

    # starting simulation loop
    for i in range( MMsys.loop_range ):

        MMsys.simmd.step( MMsys.write_frequency )

        print(i,strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        print(i,datetime.now())

        if MMsys_CPU is not None:

            state = MMsys.simmd.context.getState( getEnergy=True,getForces=True,getPositions=True )
            pos_list = state.getPositions( )
            MMsys_CPU.simmd.context.setPositions( pos_list )

            # checking if charge data needs to be collected and writing to file
            if MMsys.collect_charge_data:

                ( E , embed_charge ) = run_qmmm( QMsys , MMsys_CPU , QMsys_tare )

                fh.write( '\n' )
                fh.write( str(E) + '\t' + str(embed_charge) )
                fh.flush()

            else:

                E = run_qmmm( QMsys , MMsys_CPU , QMsys_tare )

                fh.write( '\n' )
                fh.write( str(E) )
                fh.flush()

        else:

            # checking if charge data needs to be collected and writing to file
            if MMsys.collect_charge_data:

                ( E , embed_charge ) = run_qmmm( QMsys , MMsys , QMsys_tare )

                fh.write( '\n' )
                fh.write( str(E) + '\t' + str(embed_charge) )
                fh.flush()

            else:

                E = run_qmmm( QMsys , MMsys , QMsys_tare )

                fh.write( '\n' )
                fh.write( str(E) )
                fh.flush()

        # printing energies for this iteration of the simulation
        print( str(state.getKineticEnergy()) )
        print( str(state.getPotentialEnergy()) )

        for j in range( MMsys.system.getNumForces() ):

            f = MMsys.system.getForce( j )
            print( type(f) , str(MMsys.simmd.context.getState( getEnergy=True , groups=2**j ).getPotentialEnergy()) )

    fh.close()

    # checking and creating final coordinates pdb
    if MMsys.return_system_final_pdb:
        MMsys.write_pdb( 'final' )

    print('Done!')

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

