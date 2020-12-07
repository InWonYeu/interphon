import numpy as np
from typing import List, Dict
from InterPhon.util import FilePath, File
from .unit_cell import UnitCell
from .super_cell import SuperCell
from .pre_check import PreArgument


class PreProcess(object):
    """
    Pre process class to control pre-process.
    The information required in pre-process is stored in instance variables which are determined by
    'set_user_arg', 'set_unit_cell', and 'set_super_cell' methods.
    """
    def __init__(self, user_arg=PreArgument(), unit_cell=UnitCell(), super_cell=SuperCell()):
        """
        Constructor of PreProcess class.

        :param user_arg: (instance) of PreArgument class.
        :param unit_cell: (instance) of UnitCell class.
        :param super_cell: (instance) of SuperCell class.
        """
        self.__user_arg = user_arg
        self.__unit_cell = unit_cell
        self.__super_cell = super_cell

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
        Method of PreProcess class.
        Process to set the PreArgument instance from the information given by user.

        usage:
        " >>> instance_of_PreProcess.set_user_argument(list_args=arguments)"

        :param dict_args: (Dict) Arguments given by user.
        :return: (None)
        """
        self.user_arg.set_user_argument(dict_args)
        self.user_arg.check_user_argument()

    def set_unit_cell(self, in_file: FilePath, code_name: str = 'vasp') -> None:
        """
        Method of PreProcess class.
        Process to set the variables of UnitCell instance by reading input file of DFT programs.

        usage:
        " >>> instance_of_PreProcess.set_unit_cell(in_file='./POSCAR')"

        :param in_file: (str) Path of the input file.
        :param code_name: (str) Specification of the file-format.
        :return: (None)
        """
        self.unit_cell.read_unit_cell(in_file, code_name=code_name)

    def set_super_cell(self, out_file: FilePath, comment: str = 'Supercell',
                       write_file: bool = True, code_name: str = 'vasp') -> File:
        """
        Method of PreProcess class.
        Process to set the variables of SuperCell instance from the information of UnitCell and UserArgument instances.

        usage:
        " >>> instance_of_PreProcess.set_super_cell(out_file='./SPOSCAR', comment='Supercell', write_file=True)"

        :param out_file: (str) Name of DFT input file.
        :param comment: (str) Comment of DFT input file.
        :param write_file: (bool) Write (True) the DFT input file from the information of SuperCell instance
        or not (False)
        :param code_name: (str) Specification of the file-format.
        :return: (File)
        """
        self.super_cell.set_super_cell(self.unit_cell, self.user_arg)
        if write_file is True:
            self.super_cell.write_unit_cell(out_file, comment=comment, code_name=code_name)

    def write_displace_cell(self, out_file: FilePath, code_name: str = 'vasp') -> File:
        """
        Method of PreProcess class.
        Process to write the displaced SuperCell instance into DFT input file format.

        usage:
        " >>> instance_of_PreProcess.write_displace_cell(out_file='./DPOSCAR')"

        :param out_file: (str) Name of DFT input file.
        :param code_name: (str) Specification of the file-format.
        :return: (File)
        """
        _dis_super_cell = self.super_cell
        _current_position = self.super_cell.atom_cart.copy()
        _enlarge = 1
        for ind, value in enumerate(self.user_arg.periodicity):
            if value:
                _enlarge = _enlarge * self.user_arg.enlargement[ind]

        for i, ind_T in enumerate(self.unit_cell.atom_true):
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
