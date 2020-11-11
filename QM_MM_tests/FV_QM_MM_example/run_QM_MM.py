from __future__ import print_function
import sys
# append path to QM/MM class library
sys.path.append('../../lib')
# append path to Fixed-Voltage class library
sys.path.append('/home/mcdanielgroup/data/Jesse/Ferrocene_electrode/Fixed_Voltage_OpenMM/lib')
#********* import QMclass
from QM_classes import *
#********* import additional subroutines
from QM_MM_ixn import *
#********* import MMclasses
from MM_class_base import *
from MM_classes import *
from MM_classes_FV import *
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
import argparse

# collecting and parsing user argument for input file
parser = argparse.ArgumentParser( description='input file' )
parser.add_argument( 'input_file' , type=str )
args = parser.parse_args()

# parsing arguments from input file
input_args = parse_input_file( args.input_file )

# input files for setting up OpenMM
pdb_list = [x.strip() for x in input_args.pop( 'pdb_list' , False ).split(',')]
residue_xml_list = [x.strip() for x in input_args.pop( 'residue_xml_list' , False ).split(',')]
ff_xml_list = [x.strip() for x in input_args.pop( 'ff_xml_list' , False ).split(',')]

#print( 'pdb_list' , pdb_list )
#print( 'residue_xml' , residue_xml_list )
#print( 'ff_xml_list' , ff_xml_list )
#sys.exit()

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit(2000)

# QM atoms and MM atoms with analytic Coulomb embedding
QMregion_list=(8152,8153,8154,8155,8156,8157)


#************************** download SAPT-FF force field files from github
url1 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt.xml"
url2 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt_residues.xml"
#url3 = "https://raw.github.com/jmcdaniel43/SAPT_force_field_OpenMM/master/sapt_exclusions.py"
filename1, headers = urllib.request.urlretrieve(url1, "ffdir/sapt.xml")
filename2, headers = urllib.request.urlretrieve(url2, "ffdir/sapt_residues.xml")
#filename3, headers = urllib.request.urlretrieve(url3, "sapt_exclusions.py")

# add extra CustomNonbonded force parameters to .xml file
add_CustomNonbondedForce_SAPTFF_parameters( xml_base = "ffdir/sapt.xml" , xml_param = "ffdir/graph_customnonbonded.xml" , xml_combine = "ffdir/sapt_add.xml" )


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

MMsys=MM_FixedVoltage_QMMM( pdb_list , residue_xml_list , ff_xml_list , **input_args )

# if periodic residue, call this
MMsys.set_periodic_residue(True)


#***********  Initialze OpenMM API's, this method creates simulation object
MMsys.set_platform( input_args['platform'] )

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
MMsys.generate_exclusions( water_name = 'HOH' , flag_hybrid_water_model = False , flag_SAPT_FF_exclusions = True  )

# Umbrella potential on QM atoms
#MMsys.setumbrella( 'N2', 'grph', 'C100', 2000.0 , 0.4 )   # molecule1, molecule2, atom2,  k (kJ/mol/nm^2) , r0 nm

# Fixed Voltage Electrostatics ...
MMsys.Poisson_solver_fixed_voltage( Niterations=3 )


# *********************************************************************
#                     Create QM system object
#**********************************************************************

QMsys = QM( **input_args )
QMsys_tare = None

# checking if a QMsys_tare needs to be initialized
if QMsys_tare is None and QMsys.qmmm_tare:

    tare_input_args = copy( input_args )
    tare_input_args['qmmm_ewald'] = 'False'
    QMsys_tare = QM( **tare_input_args )

#**********************************************************************
#                     QM/MM Energy
#**********************************************************************

( E , Q ) = qmmm_energy( QMsys , MMsys , QMregion_list = QMregion_list , QMsys_tare = QMsys_tare , collect_charge_data = input_args['collect_charge_data'] )

# generating output file
fh = open( input_args['out_dir'] , 'w')
fh.write( str(E) )
fh.close()

sys.exit()

