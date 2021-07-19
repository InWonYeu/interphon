import numpy as np
from InterPhon.util import FilePath, File
from InterPhon.util import Symmetry2D
from InterPhon import error
from .unit_cell import UnitCell
from .super_cell import SuperCell
from .pre_check import PreArgument


class PreProcess(object):
    """
    Pre process class to control pre-process.
    After pre-process is successfully done, the files of displaced supercells and a pre-process record file (**pre_process.yaml**) will be generated.
    The information required in pre-process is stored in the instance variables of this class.
    The instance variables are set by the :class:`core.PreProcess.set_user_arg`, :class:`core.PreProcess.set_unit_cell`, and :class:`core.PreProcess.set_super_cell` methods.

    :param user_arg: Instance of PreArgument class
    :type user_arg: :class:`core.PreArgument`
    :param unit_cell: Instance of UnitCell class
    :type unit_cell: :class:`core.UnitCell`
    :param super_cell: Instance of SuperCell class
    :type super_cell: :class:`core.SuperCell`
    """
    def __init__(self, user_arg=PreArgument(),
                 unit_cell=UnitCell(),
                 super_cell=SuperCell()):
        """
        Constructor of PreProcess class.
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

    def set_user_arg(self, dict_args: dict) -> None:
        """
        Set a **PreArgument** instance from the information given by user.

        :param dict_args: Argument dictionary given by user
        :type dict_args: dict
        """
        self.user_arg.set_user_argument(dict_args)
        self.user_arg.check_user_argument()

    def set_unit_cell(self, in_file: FilePath,
                      code_name: str = 'vasp') -> None:
        """
        Set a **UnitCell** instance by reading input file of DFT programs.

        :param in_file: Path of DFT input file to read
        :type in_file: str
        :param code_name: Indicate a DFT program corresponding to the input file, defaults to vasp
        :type code_name: str
        """
        self.unit_cell.read_unit_cell(in_file, code_name=code_name)

    def set_super_cell(self, out_file: FilePath,
                       comment: str = 'Supercell',
                       write_file: bool = True,
                       code_name: str = 'vasp') -> File:
        """
        Set a **SuperCell** instance from the information stored in its **UnitCell** and **UserArgument** instances.

        :param out_file: Name of DFT input file to write
        :type out_file: str
        :param comment: Comment to be written in the DFT input file, defaults to Supercell
        :type comment: str
        :param write_file: Write (`True`) or not (`False`) the DFT input file from the information stored in its SuperCell instance, defaults to `True`
        :type write_file: bool
        :param code_name: Specification of the file-format by a DFT program, defaults to vasp
        :type code_name: str
        """
        self.super_cell.set_super_cell(self.unit_cell, self.user_arg)
        self.super_cell.set_super_ind_true(self.unit_cell, self.user_arg)
        if write_file is True:
            self.super_cell.write_unit_cell(out_file, comment=comment, code_name=code_name)

    def write_displace_cell(self, out_file: FilePath,
                            code_name: str = 'vasp',
                            sym_flag: bool = True) -> File:
        """
        Write displaced supercell files in a format of DFT input file.

        :param out_file: Name of DFT input file to write
        :type out_file: str
        :param code_name: Specification of the file-format by a DFT program, defaults to vasp
        :type code_name: str
        :param sym_flag: Specify whether to use symmetry operation, defaults to `True`
        :type sym_flag: bool
        """
        _dis_super_cell = self.super_cell
        _dis_super_cell.selective = False
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
