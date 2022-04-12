import sys
import numpy as np
import math

#********* this file contains routines for gathering vext info from OpenMM C++ output, and computing reference vext info for benchmarking
# this is minimum image code consistent with subroutine
# ReferenceForce::getDeltaRPeriodic in ReferenceForce.cpp from OpenMM 
def r_minimum_image_pbc( r1 , r2 , box ):
    box_a = np.array( box[0] )
    box_b = np.array( box[1] )
    box_c = np.array( box[2] )

    dr = np.array( r2 ) - np.array( r1 )
    dr = dr - box_c * math.floor( dr[2]/ box[2][2] + 0.5) 
    dr = dr - box_b * math.floor( dr[1]/ box[1][1] + 0.5) 
    dr = dr - box_a * math.floor( dr[0]/ box[0][0] + 0.5) 

    rmag = math.sqrt( np.dot( dr , dr ) )
    return rmag


def compute_vext_ref( simmd , nbondedForce , pme_grid_a , pme_grid_b , pme_grid_c , QMatoms=() ):
    # get particle positions
    state = simmd.context.getState(getEnergy=False,getForces=False,getVelocities=False,getPositions=True)
    positions = state.getPositions()

    boxvec_temp = simmd.topology.getPeriodicBoxVectors()
    # get rid of units, should be nm ...
    boxvec=[]
    for i in range(3):
        temp=[]
        for j in range(3):
            temp.append( boxvec_temp[i][j]._value )
        boxvec.append( temp )

    small=1E-7
    if abs(boxvec[0][1]) > small or abs(boxvec[0][2]) > small or abs(boxvec[1][0]) > small or abs(boxvec[1][2]) > small or abs(boxvec[2][0]) > small or abs(boxvec[2][1]) > small :
        sys.exit('vext reference potential is only coded in for orthorhombic simulation box \n')
    
    # get charges from NonBondedforce 
    charges=[]
    for res in simmd.topology.residues():
        for i in range(len(res._atoms)):
            index = res._atoms[i].index
            (q, sig, eps) = nbondedForce.getParticleParameters(index)
            charges.append( q._value )


    #print( 'box' , boxvec[0][0] , boxvec[1][1] , boxvec[2][2] )
    #for i in range(len(positions)):
    #    print( i , charges[i] , positions[i][0] , positions[i][1] , positions[i][2] )


    # distance in nm, convert to Bohr, and then convert potential to kJ/mol/charge unit
    conv = 1.0 / 18.8973 * 2625.5
 
    # loop over grid points and compute contribution of Coulomb interaction from each atom
    box_a = np.array( boxvec[0] )
    box_b = np.array( boxvec[1] )
    box_c = np.array( boxvec[2] )    
 
    vext_reference =[[ [ 0.0 for col in range(pme_grid_a)] for col in range(pme_grid_b)] for row in range(pme_grid_c)]

    for i in range(pme_grid_a):
        for j in range(pme_grid_b):
            for k in range(pme_grid_c):
  
                # absolute position of grid point
                r_grid= float(i)/float(pme_grid_a)*box_a + float(j)/float(pme_grid_b)*box_b + float(k)/float(pme_grid_c)*box_c

                vext_grid=0.0
                # loop over atoms
                for i_atom in range(len(charges)):
                    if i_atom in QMatoms:
                       continue
                    else:
                        r_atom=[]
                        r_atom.append( positions[i_atom][0]._value )
                        r_atom.append( positions[i_atom][1]._value )
                        r_atom.append( positions[i_atom][2]._value )      
    
                        # pbc
                        rmag = r_minimum_image_pbc( r_atom , r_grid , boxvec )
                        vext_grid += charges[i_atom] / rmag * conv 

                vext_reference[i][j][k] = vext_grid     

    return vext_reference




def pull_vext_data( data , pme_grid_a , pme_grid_b , pme_grid_c ):
    # these are headers to search for start of grid info.
    # these are hard-coded in to C++ code ...
    header_recip="reciprocal space external potential"
    header_real ="real space plus reciprocal space external potential"

    # find reciprocal space vext first ...
    index=0
    flag=1
    for line in data:
        index+=1
        if header_recip in line:
            flag=0
            break

    if flag == 1 :
        sys.exit("couldn't find reciprocal space section in data!\n")

    vext_recip=[]
    vext_total=[]

    # read in grid indices to make sure they match what we expect
    for i in range(pme_grid_a):
        for j in range(pme_grid_b):
            for k in range(pme_grid_c): 
                string = data[index]
                temp = string.split()
                if i != int(temp[0]) or j != int(temp[1]) or k != int(temp[2]):
                    sys.exit('pme grid indices dont match in data read!')
              
                vext_recip.append( float(temp[3]) )
                index+=1

           
    # now find real space vext ...

    index=0
    flag=1
    for line in data:
        index+=1
        if header_real in line:
            flag=0
            break

    if flag == 1 :
        sys.exit("couldn't find real space section in data!\n")


    # read in grid indices to make sure they match what we expect
    for i in range(pme_grid_a):
        for j in range(pme_grid_b):
            for k in range(pme_grid_c):
                string = data[index]
                temp = string.split()
                if i != int(temp[0]) or j != int(temp[1]) or k != int(temp[2]):
                    sys.exit('pme grid indices dont match in data read!')

                vext_total.append( float(temp[3]) )
                index+=1

 

    return vext_recip , vext_total

