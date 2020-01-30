from __future__ import print_function
import io
import ctypes
import os, sys
from redirect_stdout import stdout_redirector
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#******** this is module that goes with sapt force field files to generate exclusions
from sapt_exclusions import *
#***************************
#******* this module has routines for texting vext
from routines import *

#pdb = PDBFile('spce_spherical.pdb')
pdb = PDBFile('spce_spherical_biggerbox.pdb')

temperature= 300*kelvin

# PME grid dimensions and alpha
# choice of alpha:  For n=43 grid, 60 Angstrom box, OpenMM chooses alpha= 2.389328 nm^-1
pme_grid_a = 60
pme_grid_b = 60
pme_grid_c = 60
pme_alpha = 3.0 # nm^-1

integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
integ_md.setMaxDrudeDistance(0.02)  # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)

pdb.topology.loadBondDefinitions('residues.xml')
pdb.topology.createStandardBonds();
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('spce.xml')
modeller.addExtraParticles(forcefield)

# remember, indexing starts at zero ...
QMatoms = (120, 121, 122)
# Add tuple of QM atoms to topology
modeller.topology.addQMatoms( QMatoms )

system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True)

nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
nbondedForce.setNonbondedMethod(NonbondedForce.PME)
customNonbondedForce.setNonbondedMethod(min(nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))
customNonbondedForce.setUseLongRangeCorrection(True)

nbondedForce.setPMEParameters( pme_alpha , pme_grid_a , pme_grid_b , pme_grid_c ) # alpha, nx, ny, nz

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

properties = {'ReferenceVextGrid': 'true'}
platform = Platform.getPlatformByName('Reference')

simmd = Simulation(modeller.topology, system, integ_md, platform, properties)
simmd.context.setPositions(modeller.positions)


#************************************************
#         IMPORTANT: generate exclusions for SAPT-FF
#
sapt_exclusions = sapt_generate_exclusions(simmd,system,modeller.positions)
#
#************************************************

# call Energy and Force to test vext
# redirect stdout from C++ into python buffer

#state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)

state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,getVext_grids=True, getPME_grid_positions=True)

vext_tot = state.getVext_grid()
PME_grid_positions = state.getPME_grid_positions()

#sys.exit()

# print PME grid
#print(" PME grid \n")
#for i in range(pme_grid_a):
#    for j in range(pme_grid_b):
#        for k in range(pme_grid_c):
#            index = i * pme_grid_b * pme_grid_c + j * pme_grid_c + k
#            print( i , j , k , PME_grid_positions[index][0] ,  PME_grid_positions[index][1] ,  PME_grid_positions[index][2] )


#libc = ctypes.CDLL(None)
#f = io.BytesIO()
#with stdout_redirector(f):
#    libc.puts( state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True) )

#string = "{0}".format(f.getvalue().decode('utf-8'))
#data = string.splitlines()

# get vext on grid from data string
#vext_recip , vext_tot = pull_vext_data( data , pme_grid_a , pme_grid_b , pme_grid_c )

# compute reference vext
vext_ref = compute_vext_ref( simmd , nbondedForce , pme_grid_a , pme_grid_b , pme_grid_c , QMatoms )


# print PME computed reciprocal space and real+reciprocal vext on grid points vs benchmark
print( "Evaluation of vext on PME grid: i, j , k ,  PME_recip , PME_real + PME_recip , benchmark : \n")
for i in range(pme_grid_a):
    for j in range(pme_grid_b):
        for k in range(pme_grid_c):
            index = i * pme_grid_b * pme_grid_c + j * pme_grid_c + k
            
            #print( i , j , k , vext_recip[index] , vext_tot[index] , vext_ref[i][j][k] )
            print( vext_tot[index] , vext_ref[i][j][k] )

exit()
