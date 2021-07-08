import numpy as np
from .unit_cell import UnitCell
from .pre_check import PreArgument
from InterPhon.util import MatrixLike, AtomType, SelectIndex


class SuperCell(UnitCell):
    """
    Super cell class to construct an super cell object in a standardized format.
    This child class is inherited from the :class:`core.UnitCell` parent class.
    The information about atomic structure of super cell is stored in the instance variables of this class.
    The instance variables are set by the :class:`core.SuperCell.set_super_cell` and :class:`core.SuperCell.set_super_ind_true` methods.

    :param lattice_matrix: '(3, 3) size' matrix for the lattice vectors of unit cell, defaults to None
    :type lattice_matrix: np.ndarray[float]
    :param atom_type: List of atom type in unit cell, defaults to None
    :type atom_type: AtomType
    :param num_atom: Matrix of the number of atoms of each atom type, defaults to None
    :type num_atom: np.ndarray[int]
    :param selective: Selective dynamics (True) or not (False), defaults to `False`
    :type selective: bool
    :param coordinate: 'direct' or 'cartesian' coordinate, defaults to None
    :type coordinate: str
    :param atom_cart: '(total number of atoms, 3) size' matrix for atom positions in cartesian coordinate, defaults to None
    :type atom_cart: np.ndarray[float]
    :param atom_true: Index of selected atoms which are allowed to move, defaults to None
    :type atom_true: SelectIndex
    :param xyz_true: (SelectIndex) Index of the atoms which are allowed to move for each x, y, z cartesian direction, defaults to None
    :type xyz_true: SelectIndex
    :param mass_true: '(3 * number of selected atoms, ) size' matrix for mass of selected atoms, defaults to None
    :type mass_true: np.ndarray[float]
    """
    def __init__(self, lattice_matrix: np.ndarray = None,
                 atom_type: AtomType = None,
                 num_atom: np.ndarray = None,
                 selective: bool = False,
                 coordinate: str = None,
                 atom_cart: np.ndarray = None,
                 atom_true: SelectIndex = None,
                 xyz_true: SelectIndex = None,
                 mass_true: np.ndarray = None):
        """
        Constructor of SuperCell class.
        """
        super(SuperCell, self).__init__(lattice_matrix,
                                        atom_type,
                                        num_atom,
                                        selective,
                                        coordinate,
                                        atom_cart,
                                        atom_true,
                                        xyz_true,
                                        mass_true)

    def set_super_cell(self, unit_cell: UnitCell,
                       user_arg: PreArgument) -> None:
        """
        Set the variables of **SuperCell** instance from the information stored in given **UnitCell** and **UserArgument** instances.

        :param unit_cell: Instance of UnitCell class
        :type unit_cell: :class:`core.UnitCell`
        :param user_arg: Instance of PreArgument class
        :type user_arg: :class:`core.PreArgument`
        """
        self.lattice_matrix = np.asfarray([unit_cell.lattice_matrix[i, 0:3] * user_arg.enlargement[i] for i in range(3)])

        _enlarge = 1
        for ind, value in enumerate(user_arg.periodicity):
            if value:
                _enlarge = _enlarge * user_arg.enlargement[ind]

        self.atom_type = []
        for atom_type in unit_cell.atom_type:
            self.atom_type.extend([atom_type for _ in range(_enlarge)])

        self.num_atom = unit_cell.num_atom * _enlarge
        self.coordinate = unit_cell.coordinate

        _super_direct = np.empty([_enlarge, 3])
        k = 0
        for x in range(user_arg.enlargement[0]):
            for y in range(user_arg.enlargement[1]):
                for z in range(user_arg.enlargement[2]):
                    _super_direct[k, 0:3] = np.array([float(x) / user_arg.enlargement[0],
                                                      float(y) / user_arg.enlargement[1],
                                                      float(z) / user_arg.enlargement[2]])
                    k = k + 1

        self.atom_cart = np.empty([unit_cell.atom_cart.shape[0] * _enlarge, unit_cell.atom_cart.shape[1]])
        for ind, value in enumerate(unit_cell.atom_cart):
            self.atom_cart[_enlarge * ind: _enlarge * (ind + 1), 0:3] = value + \
                                                                        np.dot(_super_direct, self.lattice_matrix)

    def set_super_ind_true(self, unit_cell: UnitCell,
                           user_arg: PreArgument) -> None:
        """
        Set the **SuperCell** index of selected atoms from the information stored in given **UnitCell** and **UserArgument** instances.

        :param unit_cell: Instance of UnitCell class
        :type unit_cell: :class:`core.UnitCell`
        :param user_arg: Instance of PreArgument class
        :type user_arg: :class:`core.PreArgument`
        """
        _enlarge = 1
        for ind, value in enumerate(user_arg.periodicity):
            if value:
                _enlarge = _enlarge * user_arg.enlargement[ind]

        self.selective = True

        self.atom_true = []
        for ind_true in unit_cell.atom_true:
            self.atom_true.extend([ind_true * _enlarge + ind for ind in range(_enlarge)])

        self.xyz_true = []
        for ind_true in self.atom_true:
            self.xyz_true.extend([ind_true for i in range(3)])
