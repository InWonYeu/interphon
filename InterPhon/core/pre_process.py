import numpy as np
from typing import List, Dict
from InterPhon.util import FilePath, File
from InterPhon.util import Symmetry2D
from InterPhon import error
from .unit_cell import UnitCell
from .super_cell import SuperCell
from .pre_check import PreArgument


class PreProcess(object):
    """
    Pre process class to control pre-process.
    After pre-process is successfully done, the files of displaced supercells and a pre-process record file (pre_process.yaml) will be generated.
    The information required in pre-process is stored in the instance variables of this class.
    The instance variables are set by the 'set_user_arg', 'set_unit_cell', and 'set_super_cell' methods.
    """
    def __init__(self, user_arg=PreArgument(),
                 unit_cell=UnitCell(),
                 super_cell=SuperCell()):
        """
        Constructor of PreProcess class.

        :param user_arg: (instance) of PreArgument class.
        :param unit_cell: (instance) of UnitCell class.
        :param super_cell: (instance) of SuperCell class.
        """
        self.__user_arg = user_arg
        self.__unit_cell = unit_cell
        self.__super_cell = super_cell
        self.sym = Symmetry2D(self.unit_cell, self.super_cell, self.user_arg)

    @property
    def user_arg(self):
        return self.__user_arg

    @user_arg.setter
    def user_arg(self, _user_arg):
        if isinstance(_user_arg, PreArgument):
            self.__user_arg = _user_arg
        else:
            ValueError(
                "'{0}' should be the instance of <class 'InterPhon.core.pre_check.PreArgument'>".format(_user_arg))

    @property
    def unit_cell(self):
        return self.__unit_cell

    @unit_cell.setter
    def unit_cell(self, _unit_cell):
        if isinstance(_unit_cell, UnitCell):
            self.__unit_cell = _unit_cell
        else:
            ValueError(
                "'{0}' should be the instance of <class 'InterPhon.core.unit_cell.UnitCell'>".format(_unit_cell))

    @property
    def super_cell(self):
        return self.__super_cell

    @super_cell.setter
    def super_cell(self, _super_cell):
        if isinstance(_super_cell, SuperCell):
            self.__super_cell = _super_cell
        else:
            ValueError(
                "'{0}' should be the instance of <class 'InterPhon.core.super_cell.SuperCell'>".format(_super_cell))

    def set_user_arg(self, dict_args: Dict) -> None:
        """
        Instance method of PreProcess class.
        Set a PreArgument instance from the information given by user.

        usage:
        " >>> instance_of_PreProcess.set_user_argument(dict_args=arguments)"

        :param dict_args: (Dict) Argument dictionary given by user.
        :return: (None)
        """
        self.user_arg.set_user_argument(dict_args)
        self.user_arg.check_user_argument()

    def set_unit_cell(self, in_file: FilePath, code_name: str = 'vasp') -> None:
        """
        Instance method of PreProcess class.
        Set a UnitCell instance by reading input file of DFT programs.

        usage:
        " >>> instance_of_PreProcess.set_unit_cell(in_file='./POSCAR', code_name='vasp')"

        :param in_file: (str) Path of DFT input file to read.
        :param code_name: (str) Indicate a DFT program corresponding to the input file.
        :return: (None)
        """
        self.unit_cell.read_unit_cell(in_file, code_name=code_name)

    def set_super_cell(self, out_file: FilePath, comment: str = 'Supercell',
                       write_file: bool = True, code_name: str = 'vasp') -> File:
        """
        Instance method of PreProcess class.
        Set a SuperCell instance from the information stored in its UnitCell and UserArgument instances.

        usage:
        " >>> instance_of_PreProcess.set_super_cell(out_file='./SPOSCAR', comment='Supercell', write_file=True, code_name='vasp')"

        :param out_file: (str) Name of DFT input file to write.
        :param comment: (str) Comment to be written in the DFT input file.
        :param write_file: (bool) Write (True) or not (False) the DFT input file from the information stored in its SuperCell instance.
        :param code_name: (str) Specification of the file-format by a DFT program.
        :return: (File)
        """
        self.super_cell.set_super_cell(self.unit_cell, self.user_arg)
        self.super_cell.set_super_ind_true(self.unit_cell, self.user_arg)
        if write_file is True:
            self.super_cell.write_unit_cell(out_file, comment=comment, code_name=code_name)

    def write_displace_cell(self, out_file: FilePath, code_name: str = 'vasp', sym_flag: bool = True) -> File:
        """
        Instance method of PreProcess class.
        Write displaced supercell files in a format of DFT input file.

        usage:
        " >>> instance_of_PreProcess.write_displace_cell(out_file='./SPOSCAR', code_name='vasp')"

        :param out_file: (str) Name of DFT input file to write.
        :param code_name: (str) Specification of the file-format by a DFT program.
        :param sym_flag: (bool) Specify whether to use symmetry operation.
        :return: (File)
        """
        _dis_super_cell = self.super_cell
        _current_position = self.super_cell.atom_cart.copy()

        _enlarge = 1
        for ind, value in enumerate(self.user_arg.periodicity):
            if value:
                _enlarge = _enlarge * self.user_arg.enlargement[ind]

        try:
            if sym_flag:
                self.sym = Symmetry2D(self.unit_cell, self.super_cell, self.user_arg)
                _, _, _ = self.sym.search_point_group()
        except error.Cannot_Search_Point_Group as e:
            print("look-up table: ", e.value)
            print(e)
            sym_flag = False

        if sym_flag:
            _, _, _, _ = self.sym.search_image_atom()
            self.sym.search_self_image_atom()
            self.sym.search_independent_displacement()
            self.sym.gen_additional_displacement()

            k = 0
            for i, ind_T in enumerate(self.sym.require_atom):
                _dis_super_cell.atom_cart = _current_position.copy()

                _displace = [self.sym.independent_by_single_displacement_cart[i][0]]
                _additional_displace = self.sym.independent_additional_displacement_cart[i]
                if _additional_displace:
                    _displace.extend(_additional_displace)

                for j, displace in enumerate(_displace):
                    # Forward Displacement
                    _dis_super_cell.atom_cart[_enlarge * self.unit_cell.atom_true[ind_T], 0:3] = \
                        _current_position[_enlarge * self.unit_cell.atom_true[ind_T],
                        0:3] + self.user_arg.displacement * displace
                    _dis_super_cell.write_unit_cell(out_file + '-{0:0>4}'.format(2 * k + 1),
                                                    comment='Forward Displacement', code_name=code_name)

                    # Backward Displacement
                    _dis_super_cell.atom_cart[_enlarge * self.unit_cell.atom_true[ind_T], 0:3] = \
                        _current_position[_enlarge * self.unit_cell.atom_true[ind_T],
                        0:3] - self.user_arg.displacement * displace
                    _dis_super_cell.write_unit_cell(out_file + '-{0:0>4}'.format(2 * k + 2),
                                                    comment='Backward Displacement', code_name=code_name)

                    k += 1

        else:
            for i, ind_T in enumerate(self.unit_cell.atom_true):
                _dis_super_cell.atom_cart = _current_position.copy()
                for j, displace in enumerate(np.eye(3, dtype=float)):

                    # Forward Displacement
                    _dis_super_cell.atom_cart[_enlarge * ind_T, 0:3] = \
                        _current_position[_enlarge * ind_T, 0:3] + self.user_arg.displacement * displace
                    _dis_super_cell.write_unit_cell(out_file + '-{0:0>4}'.format(6 * i + 2 * j + 1),
                                                    comment='Forward Displacement', code_name=code_name)

                    # Backward Displacement
                    _dis_super_cell.atom_cart[_enlarge * ind_T, 0:3] = \
                        _current_position[_enlarge * ind_T, 0:3] - self.user_arg.displacement * displace
                    _dis_super_cell.write_unit_cell(out_file + '-{0:0>4}'.format(6 * i + 2 * j + 2),
                                                    comment='Backward Displacement', code_name=code_name)

    def __deepcopy__(self, memodict: dict = {}) -> object:
        import copy

        cls = self.__class__
        result = cls.__new__(cls)
        memodict[id(self)] = result
        for key, value in self.__dict__.items():
            setattr(result, key, copy.deepcopy(value, memodict))
        return result
