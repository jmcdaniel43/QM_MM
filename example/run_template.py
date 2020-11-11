from __future__ import print_function
import sys

# append library path
sys.path.append('../lib/')

#********* import QMclass
from QM_classes import *

#********* import parent/child MMclass
from MM_class_base import *
from MM_classes import *

#********* import additional subroutines
from QM_MM_ixn import *

# other stuff
from sys import stdout
import argparse

# collecting and parsing user argument for input file
parser = argparse.ArgumentParser( description='input file' )
parser.add_argument( 'input_file' , type=str )
args = parser.parse_args()

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit( 2000 )

# parsing arguments from input file
input_args = parse_input_file( args.input_file )

# *********************************************************************
#                     Create MM system object
#**********************************************************************

# input files for setting up OpenMM
pdb_list = [ input_args.pop( 'pdb_list' , False ) , ]
residue_xml_list = [ input_args.pop( 'residue_xml_list' , False ) , ]
ff_xml_list = [ input_args.pop( 'ff_xml_list' , False ) , ]

MMsys = MM_QMMM( pdb_list , residue_xml_list , ff_xml_list , **input_args )
# set Platform
MMsys.set_platform( input_args['platform'] )
# set exclusions for SAPT-FF force field
MMsys.generate_exclusions(water_name = 'HOH', flag_hybrid_water_model = False ,  flag_SAPT_FF_exclusions = True )

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

( E , Q ) = qmmm_energy( QMsys , MMsys , QMsys_tare = QMsys_tare , collect_charge_data = input_args['collect_charge_data'] )

# generating output file
fh = open( input_args['out_dir'] , 'w')
fh.write( str(E) )
fh.close()

sys.exit()
