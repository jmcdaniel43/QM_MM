#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base class for defining systems and subsystems.
"""


class System:
    """
    The System base class for defining systems and subsystems.

    Parameters
    ----------
    residue_part_list: list of list of int, Optional, default=None
        The indices of all particles in the system grouped by residue.
    group_part_dict: dict of list of list of int, Optional, default=None
        The indices of particles in the system grouped by a given key.
        Example keys include "qm_atom", "qm_drude", and "analytic".
    element_symbols: list of str, Optional, default=None
        The element symbols of all particles in the system.
    charges: NumPy Array object, Optional, default=None
        The charges of all particles in the system, in proton charge
        units.
    positions: NumPy Array object, Optional, default=None
        The positions of all particles in the system, in Angstroms.
    box: NumPy Array object, Optional, default=None
        The box vectors defining the periodic system, in Angstroms.
    """
    nm_to_bohr = 18.89726
    angstrom_to_bohr = 1.889726
    hartree_to_kjmol = 2625.4996
    angstrom_to_nm = 0.1

    def __init__(
            self,
            residue_part_list=None,
            group_part_dict=None,
            element_symbols=None,
            charges=None,
            positions=None,
            box=None,
        ):
        self._residue_part_list = residue_part_list
        self._group_part_dict = group_part_dict
        self._element_symbols = element_symbols
        self._charges = charges
        self._positions = positions
        self._box = box

    @property
    def residue_part_list(self):
        """
        The indices of particles in the system grouped by residue.
        """
        return self._residue_part_list

    @residue_part_list.setter
    def residue_part_list(self, residue_part_list):
        self._residue_part_list = residue_part_list

    @property
    def group_part_dict(self):
        """
        The indices of particles in the system grouped by a given key.
        """
        return self._group_part_dict

    @group_part_dict.setter
    def group_part_dict(self, group_part_dict):
        self._group_part_dict = group_part_dict

    @property
    def element_symbols(self):
        """
        The element symbols of all particles in the system.
        """
        return self._element_symbols

    @element_symbols.setter
    def element_symbols(self, element_symbols):
        self._element_symbols = element_symbols

    @property
    def charges(self):
        """
        The charges of all particles in the system.
        """
        return self._charges

    @charges.setter
    def charges(self, charges):
        self._charges = charges
    
    @property
    def positions(self):
        """
        The positions of all particles in the system.
        """
        return self._positions

    @positions.setter
    def positions(self, positions):
        self._positions = positions

    @property
    def box(self):
        """
        The box vectors defining the periodic system.
        """
        return self._box

    @box.setter
    def box(self, box):
        self._box = box
