#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QM/MM, a method to perform single-point QM/MM calculations using the 
QM/MM/PME direct electrostatic QM/MM embedding method.
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
    basis_set = "aug-cc-pvdz"
    functional = "PBE"
    quadrature_spherical = 302
    quadrature_radial = 75
    qm_charge = 0
    qm_spin = 1
    n_threads = 24
    # Define MM subsystem inputs.
    pdb_list = ["./data/spce_box.pdb"]
    residue_xml_list = ["./data/spce_residues.xml"]
    ff_xml_list = ["./data/spce.xml"]
    platform = "CPU"
    # Define QM/MM system inputs.
    group_part_dict = {"qm_atom": [0,1,2]}
    embedding_cutoff = 0
    embedding_method = "analytic"
    qmmm_pme = True
    qmmm_pme_gridnumber = 250
    qmmm_pme_alpha = 5.0
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
        pdb_list,
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
    cutoff_list = [3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18]
    cutoff_list = [0]
    energy_list = []
    for cutoff in cutoff_list:
        qmmm_system.embedding_cutoff = cutoff
        energy = qmmm_system.single_point_calculation()
        energy_list.append(energy)
    for cutoff, energy in zip(cutoff_list, energy_list):
        print("!    {}    {}".format(cutoff, energy), flush=True)

if __name__ == "__main__":
    main()
