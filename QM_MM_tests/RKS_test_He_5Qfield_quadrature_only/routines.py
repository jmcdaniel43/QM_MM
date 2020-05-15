import sys
# append path to MM class library
sys.path.append('/home/mcdanielgroup/data/Jesse/Fixed_Voltage_OpenMM/lib/')
#********* import MMclass
from MM_classes import *


#************** testing *************************
# this is for testing, construct a list of all the MM atoms xyz positions and charge
# to give to Psi4 to test numerical quadrature
#
#  output :: MMenv -- list of [ [ q1, x1, y1, z1 ] , [ q2, x2, y2, z2 ] ... ] 
#            for all atoms in MM region (full MM environment to embed in QM calculation)
#***********************************************
def create_MM_env_full(  MMsys , QMatoms ):
    lengthconv = 10.0  # nm to angstrom
    MMenv=[]

    # get positions from openMM
    # OpenMM state object stores positions
    state = MMsys.simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
    # these are all the positions, we only want some...
    position_all = state.getPositions()

    # loop over all atoms in system
    for atom in MMsys.simmd.topology.atoms():
        # if not QMatom
        if atom.index in QMatoms:
            pass
        else:
            # get partial charge from force field.  Need this if this is MM atom in QMregion...
            (q_i, sig, eps) = MMsys.nbondedForce.getParticleParameters(atom.index)
            MMenv.append( [ q_i._value , position_all[atom.index][0]._value*lengthconv , position_all[atom.index][1]._value*lengthconv , position_all [atom.index][2]._value*lengthconv ] )  # q , x , y , z

    return MMenv

