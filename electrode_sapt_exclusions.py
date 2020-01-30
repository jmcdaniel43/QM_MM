#********** OpenMM Drivers
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import *

#*************************** README  **************************************
#  This module defines methods that introduce additional exclusions
#  for electrodes, and specific electrolyte molecules within the SAPT force field.
#  This module should thus be considered "part of" the force field,
#  as (non-electrode) exclusions are force-field specific
#**************************************************************************


class electrode_sapt_generate_exclusions:

    """
    Set up exclusions for electrodes and certain molecules in SAPT force field.
    Exclusions should be implemented separately for different molecule types.
    This method should figure out which methods to subsequently call
    based on the system to generate molecule-specific exclusions
    """

    def __init__(self, sim, system, positions, cathode_name , anode_name):
    
        #************************** List of molecule types that have special exclusions ******
        # now see what molecule types are present.  These name are hardcoded in, i don't see any way around this since special exclusions
        # have to be molecule specific
        self.watername = 'HOH'
        self.TFSIname  = 'Tf2N'
        
        #***********************************************
       

        self.nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
        self.customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
        self.drudeForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == DrudeForce][0]
        # will only have this for certain molecules
        self.custombond = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomBondForce][0]


        #************** Exclusions for electrodes *******************
        self.cathode_name = cathode_name ; self.anode_name = anode_name
        flag_cathode=0; flag_anode=0  # make sure we match electrodes
        cathode=[] ; anode=[]

        #******** cathode first...
        for res in sim.topology.residues():
            if res.name == self.cathode_name:
                flag_cathode=1
                for atom in res._atoms:
                    cathode.append( atom.index )

        if flag_cathode == 0:  # make sure we've matched  cathode!
            print(' Couldnt find cathode...please check setting of cathode_name ! ')
            sys.exit(0)

        # generate exclusions for cathode
        self.exclusion_Electrode_NonbondedForce(sim , system, cathode)

        #******** now anode ...
        for res in sim.topology.residues():
            if res.name == self.anode_name:
                flag_anode=1
                for atom in res._atoms:
                    anode.append( atom.index )

        if flag_anode == 0:  # make sure we've matched anode!
            print(' Couldnt find anode...please check setting of anode_name ! ')
            sys.exit(0)

        # generate exclusions for cathode
        self.exclusion_Electrode_NonbondedForce(sim , system, anode)

        

        # ************* Add Water exclusions
        for res in sim.topology.residues():
            if res.name == self.watername:
                self.generate_exclusions_water(sim,system)
                break

        # ************* Add TFSI exclusions, note might refer to this as Tf2N
        for res in sim.topology.residues():
            if res.name == self.TFSIname:
                self.generate_exclusions_TFSI(sim,system)
                break

        # now reinitialize to make sure changes are stored in context
        sim.context.reinitialize()
        sim.context.setPositions(positions)



    #**********************************
    # Adds exclusions for intra-electrode interactions of atoms
    # Call this method with a list of electrode atoms 'electrode',
    # and all long range interactions between atoms on this electrode will
    # be excluded
    #**********************************
    def exclusion_Electrode_NonbondedForce(self, sim, system, electrode):

    # first figure out which exclusions we already have (1-4 atoms and closer).  The code doesn't
    # allow the same exclusion to be added twice
        self.flagexclusions = {}
        for i in range(self.customNonbondedForce.getNumExclusions()):
            (particle1, particle2) = self.customNonbondedForce.getExclusionParticles(i)
            string1=str(particle1)+"_"+str(particle2)
            string2=str(particle2)+"_"+str(particle1)
            self.flagexclusions[string1]=1
            self.flagexclusions[string2]=1

    # now add exclusions for every atom pair in electrode if we don't already have them
        for i in range(len(electrode)):
            indexi = electrode[i]
            for j in range(i+1,len(electrode)):
                indexj = electrode[j]
                string1=str(indexi)+"_"+str(indexj)
                string2=str(indexj)+"_"+str(indexi)
                if string1 in self.flagexclusions and string2 in self.flagexclusions:
                    continue
                else:
                    self.customNonbondedForce.addExclusion(indexi,indexj)
                    self.nbondedForce.addException(indexi,indexj,0,1,0,True)
 

    def generate_exclusions_water(self,sim,system):
        """
        Here we are simulating a hybrid SWM4-NDP/SAPT-FF model for water
        We use SWM4-NDP for water-water interactions, and SAPT-FF for
        water-other interactions.
   
        """

        # create Interaction Groups for hybrid water model
        self.water=set()
        self.notwater=set()
        print('Creating Interaction Groups for CustomNonBonded.  These interactions will be computed between water-other, not water-water.')
        for res in sim.topology.residues():
            if res.name == self.watername:
                for i in range(len(res._atoms)-1):
                    self.water.update([res._atoms[i].index])
            else:
                for i in range(len(res._atoms)-1):
                    self.notwater.update([res._atoms[i].index])

        self.customNonbondedForce.addInteractionGroup(self.water, self.notwater)
        self.customNonbondedForce.addInteractionGroup(self.notwater, self.notwater)


    def generate_exclusions_TFSI(self,sim,system):
        """
        This creates exclusions for TFSI nonbonded interactions, and update
        Screened Drude interactions.  1-5 non-Coulomb interaction are accounted for
        using CustomBondForce
        """
        print('Creating Exclusions for TFSI')

        # map from global particle index to drudeforce object index
        particleMap = {}
        for i in range(self.drudeForce.getNumParticles()):
            particleMap[self.drudeForce.getParticleParameters(i)[0]] = i

        # can't add duplicate ScreenedPairs, so store what we already have
        flagexceptions = {}
        for i in range(self.nbondedForce.getNumExceptions()):
            (particle1, particle2, charge, sigma, epsilon) = self.nbondedForce.getExceptionParameters(i)
            string1=str(particle1)+"_"+str(particle2)
            string2=str(particle2)+"_"+str(particle1)
            flagexceptions[string1]=1
            flagexceptions[string2]=1

        # can't add duplicate customNonbonded exclusions, so store what we already have
        flagexclusions = {}
        for i in range(self.customNonbondedForce.getNumExclusions()):
            (particle1, particle2) = self.customNonbondedForce.getExclusionParticles(i)
            string1=str(particle1)+"_"+str(particle2)
            string2=str(particle2)+"_"+str(particle1)
            flagexclusions[string1]=1
            flagexclusions[string2]=1

        # add exclusions for all atom pairs on TFSI residues, and when a drude pair is
        # excluded add a corresponding screened thole interaction in its place
        for res in sim.topology.residues():
            if res.name == self.TFSIname:
                for i in range(len(res._atoms)-1):
                    for j in range(i+1,len(res._atoms)):
                        (indi,indj) = (res._atoms[i].index, res._atoms[j].index)
                        # here it doesn't matter if we already have this, since we pass the "True" flag
                        self.nbondedForce.addException(indi,indj,0,1,0,True)
                        # make sure we don't already exlude this customnonbond
                        string1=str(indi)+"_"+str(indj)
                        string2=str(indj)+"_"+str(indi)
                        if string1 in flagexclusions and string2 in flagexclusions:
                            continue
                        else:
                            self.customNonbondedForce.addExclusion(indi,indj)
                        # add thole if we're excluding two drudes
                        if indi in particleMap and indj in particleMap:
                            # make sure we don't already have this screened pair
                            if string1 in flagexceptions or string2 in flagexceptions:
                                continue
                            else:
                                drudei = particleMap[indi]
                                drudej = particleMap[indj]
                                self.drudeForce.addScreenedPair(drudei, drudej, 2.0)




