import numpy as np
from .unit_cell import UnitCell
from .pre_check import PreArgument
from InterPhon.util import MatrixLike, AtomType, SelectIndex


class SuperCell(UnitCell):
    """
    Super cell class to construct super cell object in a unified single format.
    This child class is inherited from the UnitCell parent class.
    The instance variables (information about super cell) are determined by
    'set_super_cell' and 'set_super_ind_true' methods that require UnitCell and UserArgument instances.
    """
    def __init__(self, lattice_matrix: np.ndarray = None,
                 atom_type: AtomType = None,
                 num_atom: np.ndarray = None,
                 coordinate: str = None,
                 atom_cart: np.ndarray = None,
                 atom_true: SelectIndex = None,
                 xyz_true: SelectIndex = None,
                 mass_true: np.ndarray = None):
        """
        Constructor of SuperCell class.

        :param lattice_matrix: (np.ndarray[float]) '(3, 3) size' matrix for the lattice of unit cell.
        :param atom_type: (AtomType) List of atom type in unit cell.
        :param num_atom: (np.ndarray[int]) matrix of the number of atoms of each atom type.
        :param coordinate: (str) 'direct' or 'cartesian' coordinate.
        :param atom_cart: (np.ndarray[float]) '(total number of atoms, 3) size'
                                        matrix for atom positions in cartesian coordinate.
        :param atom_true: (SelectIndex) Index of selected atoms which are allowed to move.
        :param xyz_true: (SelectIndex) Index of the atoms which are allowed to move for each x,y,z direction.
        :param mass_true: (np.ndarray[float]) '(3 * number of selected atoms, ) size' matrix for mass.
        """
        super(SuperCell, self).__init__(lattice_matrix,
                                        atom_type,
                                        num_atom,
                                        coordinate,
                                        atom_cart,
                                        atom_true,
                                        xyz_true,
                                        mass_true)

    def set_super_cell(self, unit_cell: UnitCell, user_arg: PreArgument) -> None:
        """
        Method of SuperCell class.
        Set the variables of SuperCell instance from the information of UnitCell and UserArgument instances.

        usage:
        " >>> instance_of_SuperCell.set_super_cell(unit_cell=instance_of_UnitCell, user_arg=instance_of_UserArgument)"

        :param unit_cell: (instance) of UnitCell class.
        :param user_arg: (instance) of UserArgument class.
        :return: (None)
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

    def set_super_ind_true(self, unit_cell: UnitCell, user_arg: PreArgument) -> None:
        """
        Method of SuperCell class.
        Set the SuperCell index of selected atoms in UnitCell instance
        from the information of UnitCell and UserArgument instances.

        usage:
        " >>> instance_of_SuperCell.set_super_ind_true(unit_cell=instance_of_UnitCell, user_arg=instance_of_UserArgument)"

        :param unit_cell: (instance) of UnitCell class.
        :param user_arg: (instance) of UserArgument class.
        :return: (None)
        """
        _enlarge = 1
        for ind, value in enumerate(user_arg.periodicity):
            if value:
                _enlarge = _enlarge * user_arg.enlargement[ind]

        self.atom_true = []
        for ind_true in unit_cell.atom_true:
            self.atom_true.extend([ind_true * _enlarge + ind for ind in range(_enlarge)])

        self.xyz_true = []
        for ind_true in self.atom_true:
            self.xyz_true.extend([ind_true for i in range(3)])
