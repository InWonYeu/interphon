import numpy as np
import os.path
from typing import List


def read_input_lines(structure_file: str) -> tuple:
    """
    Parser function to read VASP input file.

    :param structure_file: (str) Path of VASP input file.
    :return: (tuple) A set of lattice_matrix, atom_type, coordinate, atom_cart, atom_true, xyz_true.
    """
    try:
        with open(structure_file, 'r') as infile:
            _lines = infile.readlines()
    except IOError:
        print("\nFail to open '{0}' file".format(structure_file))
        print("Please check the path of VASP input file\n")
        raise

    # comment = _lines[0]
    lattice_ratio = float(_lines[1].split()[0])
    lattice_matrix = np.asfarray([line.split()[0:3] for line in _lines[2:5]]) * lattice_ratio

    num_atom = [int(num) for num in _lines[6].split()]
    num_total = 0
    for num in num_atom:
        num_total += num

    __atom_type = []
    for ind_atom, num in enumerate(num_atom):
        __atom_type.extend([_lines[5].split()[ind_atom] for _ in range(num)])

    coordinate = 'cartesian'
    if _lines[7].split()[0][0] in ('S', 's'):
        pos_atom = _lines[9:9 + num_total]

        if _lines[8].split()[0][0] in ('D', 'd'):
            coordinate = 'direct'
            atom = [pos.split()[0:3] for pos in pos_atom]
            atom_direct = np.asfarray(atom)
            __atom_cart = np.dot(atom_direct, lattice_matrix)
        elif _lines[8].split()[0][0] in ('C', 'c', 'K', 'k'):
            atom = [pos.split()[0:3] for pos in pos_atom]
            __atom_cart = np.asfarray(atom)

        __atom_true = [i for i, v in enumerate(pos_atom) if 'T' in v.split()]

    else:
        pos_atom = _lines[8:8 + num_total]

        if _lines[7].split()[0][0] in ('D', 'd'):
            coordinate = 'direct'
            atom = [pos.split()[0:3] for pos in pos_atom]
            atom_direct = np.asfarray(atom)
            __atom_cart = np.dot(atom_direct, lattice_matrix)
        elif _lines[7].split()[0][0] in ('C', 'c', 'K', 'k'):
            atom = [pos.split()[0:3] for pos in pos_atom]
            __atom_cart = np.asfarray(atom)

        __atom_true = [i for i, _ in enumerate(pos_atom)]  # set all the atoms as True if not selective dynamics

    # # Rearrangement of the atomic position sequence according to the atomic type
    # set_atom_type = []
    # for atom in __atom_type:
    #     if atom not in set_atom_type:
    #         set_atom_type.append(atom)
    #
    # tmp_atom = []
    # for atom in __atom_type:
    #     for ind_set_atom, set_atom in enumerate(set_atom_type):
    #         if atom == set_atom:
    #             tmp_atom.append(ind_set_atom)
    # set_atom_arg = np.argsort(np.array(tmp_atom))
    #
    # _atom_cart = __atom_cart[set_atom_arg]
    # atom_type = []
    # _atom_true = []
    # for ind_atom_arg, atom_arg in enumerate(set_atom_arg):
    #     atom_type.append(__atom_type[atom_arg])
    #     if atom_arg in __atom_true:
    #         _atom_true.append(ind_atom_arg)
    #
    # # Rearrangement of the atomic position sequence according to the z-position in ascending order
    # num_atom = [atom_type.count(atom) for atom in set_atom_type]
    # num1 = 0
    # set_atom_arg = np.array([], dtype=int)
    # for num in num_atom:
    #     set_atom_arg = np.concatenate((set_atom_arg, np.argsort(_atom_cart[num1:num1 + num, 2]) + num1), axis=0)
    #     num1 += num
    #
    # atom_cart = _atom_cart[set_atom_arg]
    # atom_true = []
    # for ind_atom_arg, atom_arg in enumerate(set_atom_arg):
    #     if atom_arg in _atom_true:
    #         atom_true.append(ind_atom_arg)

    # Define the index of xyz selective dynamics True
    xyz_true = []
    for ind_T in __atom_true:
        xyz_true.extend([ind_T for i in range(3)])

    return lattice_matrix, __atom_type, num_atom, coordinate, __atom_cart, __atom_true, xyz_true


def write_input_lines(unit_cell, comment: str) -> List[str]:
    """
    Parser function to write VASP input file.

    :param unit_cell: (instance) of UnitCell class
    :param comment: (str) Comment of VASP input file.
    :return: (List[str]) List of each line of VASP input file.
    """
    lines = ["%s" % comment + '\n', "%s" % '1.00000000000000' + '\n']

    for v in unit_cell.lattice_matrix:
        _line = " {0:>20.16f}  {1:>20.16f}  {2:20.16f}".format(v[0], v[1], v[2])
        lines.append(_line + '\n')

    ind_atom = 0
    set_atom_type = [unit_cell.atom_type[ind_atom]]
    for num in unit_cell.num_atom[:-1]:
        ind_atom += num
        set_atom_type.append(unit_cell.atom_type[ind_atom])
    lines.append("%s" % "    ".join(set_atom_type) + '\n')

    num_atom = [str(num) for num in unit_cell.num_atom]
    lines.append("%s" % "    ".join(num_atom) + '\n')

    lines.append("Cartesian" + '\n')

    for v in unit_cell.atom_cart:
        _line = " {0:>20.16f}  {1:>20.16f}  {2:20.16f}".format(v[0], v[1], v[2])
        lines.append(_line + '\n')

    return lines


def read_output_lines(force_file: str, num_super_atom: int) -> np.ndarray:
    """
    Parser function to read VASP output file in which the atomic forces are written.

    :param force_file: (str) Path of VASP output file.
    :param num_super_atom: (int) The number of atoms in super cell.
    :return: _force_matrix
    """
    try:
        with open(force_file, 'r') as infile:
            _lines = infile.readlines()
    except IOError:
        print("\nFail to open '{0}' file".format(force_file))
        print("Please check the path of VASP output file\n")
        raise

    _, _filename = os.path.split(force_file)
    _unit_convert = (1.602 * 10 ** (-19)) / 10 ** (-10)  # (eV/Angst) to (J/m)
    _force_index = None

    if _filename == 'vasprun.xml':
        _tag_atomic_force = 'forces'

        for _ind_line, _line in enumerate(_lines):
            if _tag_atomic_force in _line:
                _force_index = _ind_line
                break

        if _force_index is None:
            print("'Forces acting on atoms' is not written in '{0}'".format(force_file))
            print("Corresponding DFT calculation may be incompletely stopped")
            assert False

        _atomic_forces = _lines[_force_index + 1: _force_index + 1 + num_super_atom]
        _force_matrix = np.asfarray([atomic_force.split()[1:4] for atomic_force in _atomic_forces]) * _unit_convert

    elif _filename == 'OUTCAR':
        _tag_atomic_force = 'TOTAL-FORCE'

        for _ind_line, _line in enumerate(_lines):
            if _tag_atomic_force in _line:
                _force_index = _ind_line
                break

        if _force_index is None:
            print("'Forces acting on atoms' is not written in '{0}'".format(force_file))
            print("Check: corresponding DFT calculation must have been incompletely stopped")
            assert False

        _atomic_forces = _lines[_force_index + 2: _force_index + 2 + num_super_atom]
        _force_matrix = np.asfarray([atomic_force.split()[3:6] for atomic_force in _atomic_forces]) * _unit_convert

    return _force_matrix
