import numpy as np
from typing import List
from InterPhon.util import get_atomic_weight


def read_input_lines(structure_file: str) -> tuple:
    """
    Parser function to read Quantum ESPRESSO input file.

    :param structure_file: (str) Path of Quantum ESPRESSO input file.
    :return: (tuple) A set of lattice_matrix, atom_type, coordinate, atom_cart, atom_true, xyz_true.
    """
    try:
        with open(structure_file, 'r') as infile:
            _lines = infile.readlines()
    except IOError:
        print("\nFail to open '{0}' file".format(structure_file))
        print("Please check the path of Quantum ESPRESSO input file\n")
        raise

    num_total = None
    _lattice_index = None
    _pos_atom_index = None

    for _ind_line, _line in enumerate(_lines):
        if 'nat' in _line:
            num_total = int(_line.split()[2])
        elif 'CELL_PARAMETERS' in _line:
            _lattice_index = _ind_line
        elif 'ATOMIC_POSITIONS' in _line:
            _pos_atom_index = _ind_line

    if num_total is None:
        print("'The number of atom (tag: nat)' is not written in '{0}'".format(structure_file))
        assert False
    elif _lattice_index is None:
        print("'CELL_PARAMETERS' is not written in '{0}'".format(structure_file))
        assert False
    elif _pos_atom_index is None:
        print("'ATOMIC_POSITIONS' is not written in '{0}'".format(structure_file))
        assert False

    for _ind_line, _line in enumerate(_lines[_lattice_index + 1:]):
        if _line == '\n':
            _lattice_index += 1
        else:
            break

    for _ind_line, _line in enumerate(_lines[_pos_atom_index + 1:]):
        if _line == '\n':
            _pos_atom_index += 1
        else:
            break

    lattice_matrix = np.asfarray([_line.split()[0:3] for _line in _lines[_lattice_index + 1:_lattice_index + 4]])

    __atom_type = []
    __atom_cart = []
    __atom_true = []
    for _ind_line, _line in enumerate(_lines[_pos_atom_index + 1:_pos_atom_index + 1 + num_total]):
        __atom_type.append(_line.split()[0])
        __atom_cart.append(_line.split()[1:4])
        if '0' in _line.split():
            pass
        else:
            __atom_true.append(_ind_line)

    set_atom_type = []
    num_atom = []
    pre_atom = None
    num = 0
    for atom in __atom_type:
        if atom != pre_atom:
            pre_atom = atom
            if num != 0:
                num_atom.append(num)
            set_atom_type.append(atom)
            num = 1
        else:
            num += 1
    num_atom.append(num)

    coordinate = 'cartesian'
    __atom_cart = np.asfarray(__atom_cart)

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
    Parser function to write Quantum ESPRESSO input file.

    :param unit_cell: (instance) of UnitCell class
    :param comment: (str) Comment of Quantum ESPRESSO input file.
    :return: (List[str]) List of each line of Quantum ESPRESSO input file.
    """
    lines = ["&CONTROL" + '\n' + "/\n"]
    _line = "&SYSTEM" + '\n'
    _line += "    {0:<20}".format('ntyp') + " = {0:<}".format(len(set(unit_cell.atom_type))) + '\n'
    _line += "    {0:<20}".format('nat') + " = {0:<}".format(len(unit_cell.atom_type)) + '\n'
    _line += "    {0:<20}".format('ibrav') + " = {0:<}".format(0) + '\n'
    _line += "/\n"
    lines.append(_line)
    lines.append("&ELECTRONS" + '\n' + "/\n")
    lines.append("&IONS" + '\n' + "/\n")
    lines.append("&CELL" + '\n' + "/\n")
    lines.append('\n')

    # set_atom_type = []
    # for atom in unit_cell.atom_type:
    #     if atom not in set_atom_type:
    #         set_atom_type.append(atom)
    _line = "ATOMIC_SPECIES" + '\n'
    for atom in set(unit_cell.atom_type):
        _line += atom + " {0:<}".format(get_atomic_weight(atom)) + " {0:<}_dummy.URF".format(atom) + '\n'
    lines.append(_line + '\n')

    lines.append("K_POINTS gamma" + '\n' + '\n')

    _line = "CELL_PARAMETERS angstrom" + '\n'
    for v in unit_cell.lattice_matrix:
        _line += " {0:>20.16f}  {1:>20.16f}  {2:20.16f}".format(v[0], v[1], v[2]) + '\n'
    lines.append(_line + '\n')

    _line = "ATOMIC_POSITIONS angstrom" + '\n'
    for atom, pos_atom in zip(unit_cell.atom_type, unit_cell.atom_cart):
        _line += "{0} {1:>20.16f}  {2:>20.16f}  {3:20.16f}".format(atom, pos_atom[0], pos_atom[1], pos_atom[2]) + '\n'
    lines.append(_line + '\n')

    return lines


def read_output_lines(force_file: str, num_super_atom: int) -> np.ndarray:
    """
    Parser function to read Quantum ESPRESSO output file in which the atomic forces are written.

    :param force_file: (str) Path of Quantum ESPRESSO output file.
    :param num_super_atom: (int) The number of atoms in super cell.
    :return: _force_matrix
    """
    try:
        with open(force_file, 'r') as infile:
            _lines = infile.readlines()
    except IOError:
        print("\nFail to open '{0}' file".format(force_file))
        print("Please check the path of Quantum ESPRESSO output file\n")
        raise

    _unit_convert = (13.605662 * 1.602 * 10 ** (-19)) / (0.529177 * 10 ** (-10))  # (Ry/au) to (J/m)
    _force_index = None
    _tag_atomic_force = 'Forces acting on atoms'

    for _ind_line, _line in enumerate(_lines):
        if _tag_atomic_force in _line:
            _force_index = _ind_line
            break

    if _force_index is None:
        print("'Forces acting on atoms' is not written in '{0}'".format(force_file))
        print("Corresponding DFT calculation may be incompletely stopped")
        assert False

    _atomic_forces = _lines[_force_index + 2: _force_index + 2 + num_super_atom]
    _force_matrix = np.asfarray([atomic_force.split()[6:9] for atomic_force in _atomic_forces]) * _unit_convert

    return _force_matrix
