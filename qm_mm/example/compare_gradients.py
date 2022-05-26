#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QM/MM, a method to perform single-point QM/MM calculations using the 
QM/MM/PME direct electrostatic QM/MM embedding method.

Imports
-------
os: Standard
sys: Standard
"""
import os
import sys

sys.path.append("../../")
sys.setrecursionlimit(2000)

from qm_mm import *


def main():
    """
    """
    # Define QM subsystem inputs.
    basis_set = "def2-SVP"
    functional = "PBE"
    quadrature_spherical = 146
    quadrature_radial = 75
    qm_charge = 0
    qm_spin = 1
    n_threads = 24
    # Define MM subsystem inputs.
    pdb_lists = [
        "./data/grad_i.pdb",
        "./data/grad_fx.pdb",
        "./data/grad_fy.pdb",
        "./data/grad_fz.pdb",
    ]
    residue_xml_list = ["./data/spce_residues.xml"]
    ff_xml_list = ["./data/spce.xml"]
    platform = "CPU"
    # Define QM/MM system inputs.
    group_part_dict = {"qm_atom": [0,1,2]}
    embedding_cutoff = 9
    embedding_method = "analytic"
    qmmm_pme = True
    qmmm_pme_gridnumber = 250
    qmmm_pme_alpha = 5.0
    for i, pdb in enumerate(pdb_lists):
        # Instantiate QM subsystem.
        qm_subsystem = QMSubsystem(
            basis_set,
            functional,
            quadrature_spherical,
            quadrature_radial,
            qm_charge,
            qm_spin,
            qmmm_pme=qmmm_pme,
            n_threads=n_threads,
        )
        # Instantiate MM subsystem.
        mm_subsystem = MMSubsystem(
            [pdb],
            residue_xml_list,
            ff_xml_list,
            platform,
            qmmm_pme = qmmm_pme,
            qmmm_pme_gridnumber = qmmm_pme_gridnumber,
            qmmm_pme_alpha = qmmm_pme_alpha,
        )
        # Instantiate QM/MM system.
        qmmm_system = QMMMSystem(
            qm_subsystem,
            mm_subsystem,
            group_part_dict,
            embedding_cutoff=embedding_cutoff,
            embedding_method=embedding_method,
            qmmm_pme=qmmm_pme,
            qmmm_pme_gridnumber=qmmm_pme_gridnumber,
            qmmm_pme_alpha=qmmm_pme_alpha,
        )
        # Perform calculations.
        if i == 0:
            (energy_i, forces_i) = qmmm_system.single_point_calculation(forces=True)
            print(f"CALC FOR PME={qmmm_pme_gridnumber} SPH={quadrature_spherical} RAD={quadrature_radial}")
            print("ANALYTICAL FORCES:")
            for force in forces_i:
                print(f"{force[0]} {force[1]} {force[2]}")
            print(" ")
            forces = []
        else:
            energy = qmmm_system.single_point_calculation()
            forces.append((energy_i - energy) / 0.001 / 1.8897261)
    print(f"CALC FOR PME={qmmm_pme_gridnumber} SPH={quadrature_spherical} RAD={quadrature_radial}")
    print("NUMERICAL FORCES:")
    print(f"{forces[0]} {forces[1]} {forces[2]}")
    print(" ")
    print(f"RMSD: {(sum([(y-ya)**2 for y,ya in zip(forces, forces_i[0])])/3)**0.5}")

if __name__ == "__main__":
    main()
