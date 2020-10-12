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

#********* import additional subroutines
from QM_MM_ixn import *

# other stuff
from sys import stdout
import argparse

# collecting and parsing user argument for input file
parser = argparse.ArgumentParser( description='input file directory' )
parser.add_argument( 'input_dir' , type=str )
args = parser.parse_args()

# for electrode sheets, need to up recursion limit for residue atom matching...
sys.setrecursionlimit( 2000 )

# parsing arguments from input file
input_args = parse_input_file( args.input_dir )

# *********************************************************************
#                     Create MM system object
#**********************************************************************

MMsys = MM( **input_args )

# *********************************************************************
#                     Create QM system object
#**********************************************************************

QMsys = QM( **input_args )

#**********************************************************************
#                     QM/MM Simulation
#**********************************************************************

#E = run_qmmm( QMsys , MMsys )
run_simulation( QMsys , MMsys )

# generating output file
#fh = open( input_args['out_dir'] , 'w')
#fh.write( str(E) )
#fh.close()

sys.exit()
