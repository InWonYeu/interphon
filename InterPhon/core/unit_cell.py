import numpy as np
from InterPhon.util import get_atomic_weight
from InterPhon.util import MatrixLike, AtomType, SelectIndex, FilePath, File
from InterPhon.inout import vasp, aims, espresso


def to_numpy(data):
    if isinstance(data, np.ndarray):
        pass
    elif isinstance(data, list):
        data = np.array(data)
    elif isinstance(data, tuple):
        data = np.array(data)
    elif isinstance(data, str):
        data = np.array(data.strip().split())
    else:
        raise ValueError("'{0}' cannot be set to numpy array".format(data))
    return data


def to_list(data):
    if isinstance(data, np.ndarray):
        data = list(data)
    elif isinstance(data, list):
        pass
    elif isinstance(data, tuple):
        data = list(data)
    elif isinstance(data, str):
        data = data.strip().split()
    else:
        raise ValueError("'{0}' cannot be set to list".format(data))
    return data


def to_int_list(data):
    try:
        if isinstance(data, int):
            data = [data, ]
        elif isinstance(data, str):
            data = [int(val) for val in data.strip().split()]
        else:
            data = [int(val) for val in data]
    except ValueError:
        raise ValueError("The items of '{0}' cannot be converted to int".format(data))
    return data


class UnitCell(object):
    """
    Unit cell class to construct an unit cell object in a standardized format.
    The information about atomic structure of interest is stored in the instance variables of this class.
    This class interacts with input files of different DFT programs by employing the modules in :class:`inout` sub-package.

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
        Constructor of UnitCell class.
        """
        self.__lattice_matrix = lattice_matrix
        self.__atom_type = atom_type
        self.__num_atom = num_atom
        self.selective = selective
        self.__coordinate = coordinate
        self.__atom_cart = atom_cart
        self.__atom_true = atom_true
        self.__xyz_true = xyz_true
        self.__mass_true = mass_true

    @property
    def lattice_matrix(self):
        return self.__lattice_matrix

    @lattice_matrix.setter
    def lattice_matrix(self, _lattice_matrix):
        _lattice_matrix = to_numpy(_lattice_matrix)

        if _lattice_matrix.shape == (3, 3):
            self.__lattice_matrix = _lattice_matrix
        else:
            raise ValueError("The size of lattice matrix should be (3, 3)")

    @property
    def atom_type(self):
        return self.__atom_type

    @atom_type.setter
    def atom_type(self, _atom_type):
        _atom_type = to_list(_atom_type)

        if self.atom_cart is None or len(_atom_type) == self.atom_cart.shape[0]:
            self.__atom_type = _atom_type
        else:
            raise ValueError("The length of atom_type should be the total number of atoms")

    @property
    def num_atom(self):
        return self.__num_atom

    @num_atom.setter
    def num_atom(self, _num_atom):
        _num_atom = to_numpy(_num_atom)

        if self.atom_cart is None or _num_atom.sum(axis=0) == self.atom_cart.shape[0]:
            self.__num_atom = _num_atom
        else:
            raise ValueError("The sum of atoms should be the total number of atoms")

    @property
    def coordinate(self):
        return self.__coordinate

    @coordinate.setter
    def coordinate(self, _coordinate):
        _coordinate = str(_coordinate)

        if _coordinate[0] in ('D', 'd'):
            self.__coordinate = 'direct'
        elif _coordinate[0] in ('C', 'c', 'K', 'k'):
            self.__coordinate = 'cartesian'
        else:
            raise ValueError("The first character should be in ('D', 'd') for direct (fractional) coordinates "
                             "and in ('C', 'c', 'K', 'k') for cartesian coordinates")

    @property
    def atom_cart(self):
        return self.__atom_cart

    @atom_cart.setter
    def atom_cart(self, _atom_cart):
        _atom_cart = to_numpy(_atom_cart)

        if self.atom_type is None or _atom_cart.shape == (len(self.atom_type), 3):
            self.__atom_cart = _atom_cart
        else:
            raise ValueError("The size of lattice matrix should be (total number of atoms, 3)")

    @property
    def atom_direct(self):
        return np.dot(self.atom_cart, np.linalg.inv(self.lattice_matrix))

    @atom_direct.setter
    def atom_direct(self, _atom_direct):
        _atom_direct = to_numpy(_atom_direct)

        if self.atom_type is None or _atom_direct.shape == (len(self.atom_type), 3):
            self.__atom_cart = np.dot(_atom_direct, self.lattice_matrix)
        else:
            raise ValueError("The size of lattice matrix should be (total number of atoms, 3)")

    @property
    def atom_true(self):
        return self.__atom_true

    @atom_true.setter
    def atom_true(self, _atom_true):
        _atom_true = to_int_list(_atom_true)

        if self.atom_type is None:
            self.__atom_true = sorted(_atom_true)
        else:
            mask = np.isin(np.array(_atom_true), range(len(self.atom_type)))
            if np.all(mask):
                self.__atom_true = sorted(_atom_true)
            else:
                raise ValueError("Index of atoms should be in the range(0, total number of atoms)")

        # for atom_index in _atom_true:
        #     if atom_index not in range(len(self.atom_type)):
        #         raise ValueError("Index of atoms should be in the range(0, total number of atoms)")

    @property
    def xyz_true(self):
        return self.__xyz_true

    @xyz_true.setter
    def xyz_true(self, _xyz_true):
        _xyz_true = to_int_list(_xyz_true)

        if self.atom_true is None:
            self.__xyz_true = sorted(_xyz_true)
        else:
            mask = np.isin(np.array(_xyz_true), self.atom_true)
            if np.all(mask):
                self.__xyz_true = sorted(_xyz_true)
            else:
                raise ValueError("Index of atoms should be in the range(0, total number of atoms)")

    @property
    def mass_true(self):
        return self.__mass_true

    @mass_true.setter
    def mass_true(self, _mass_true):
        _mass_true = to_numpy(_mass_true)

        if self.atom_true is None:
            self.__mass_true = _mass_true
        else:
            if _mass_true.shape == (3 * len(self.atom_true),):
                self.__mass_true = _mass_true
            else:
                raise ValueError("The size of lattice matrix should be "
                                 "(3 * number of selected atoms, )")

    def initialization(self):
        """
        Initialize the instance variables.
        """
        self.__lattice_matrix = None
        self.__atom_type = None
        self.__num_atom = None
        self.selective = False
        self.__coordinate = None
        self.__atom_cart = None
        self.__atom_true = None
        self.__xyz_true = None
        self.__mass_true = None

    def read_unit_cell(self, in_file: FilePath,
                       code_name: str = 'vasp') -> None:
        """
        Set the variables of UnitCell instance by reading input file of DFT programs.

        :param in_file: Path of DFT input file to read
        :type in_file: str
        :param code_name: Indicate a DFT program corresponding to the input file, defaults to vasp
        :type code_name: str
        """
        self.initialization()
        if code_name == 'vasp':
            self.lattice_matrix, \
            self.atom_type, \
            self.num_atom, \
            self.selective, \
            self.coordinate, \
            self.atom_cart, \
            self.atom_true, \
            self.xyz_true = vasp.read_input_lines(in_file)

        elif code_name == 'espresso':
            self.lattice_matrix, \
            self.atom_type, \
            self.num_atom, \
            self.selective, \
            self.coordinate, \
            self.atom_cart, \
            self.atom_true, \
            self.xyz_true = espresso.read_input_lines(in_file)

        elif code_name == 'aims':
            self.lattice_matrix, \
            self.atom_type, \
            self.num_atom, \
            self.selective, \
            self.coordinate, \
            self.atom_cart, \
            self.atom_true, \
            self.xyz_true = aims.read_input_lines(in_file)

    def write_unit_cell(self, out_file: FilePath,
                        comment: str = 'Unknown',
                        code_name: str = 'vasp') -> File:
        """
        Write the input file of DFT programs from the information stored in an **UnitCell** instance.

        :param out_file: Name of DFT input file to write
        :type out_file: str
        :param comment: Comment to be written in the DFT input file, defaults to Unknown
        :type comment: str
        :param code_name: Specification of the file-format by a DFT program, defaults to vasp
        :type code_name: str
        """
        if code_name == 'vasp':
            _lines = vasp.write_input_lines(self, comment)
            with open(out_file, 'w') as outfile:
                outfile.write("%s" % "".join(_lines))

        elif code_name == 'espresso':
            _lines = espresso.write_input_lines(self, comment)
            with open(out_file, 'w') as outfile:
                outfile.write("%s" % "".join(_lines))

        elif code_name == 'aims':
            _lines = aims.write_input_lines(self, comment)
            with open(out_file, 'w') as outfile:
                outfile.write("%s" % "".join(_lines))

    def set_mass_true(self) -> None:
        """
        Set the instance variable (**self.mass_true**) for selected atoms which are allowed to move
        by employing the :class:`get_atomic_weight` fuction in :class:`util` sub-package.
        """
        _mass_atom = np.asfarray([get_atomic_weight(v) for v in self.atom_type]) / (6.022 * 10 ** 23) / 1000
        self.mass_true = np.asfarray(_mass_atom)[self.xyz_true]

    # a = [1, 2, 3]
    # memo = {}
    # b = copy.deepcopy(a, memo)
    # # now memo = {139907464678864: [1, 2, 3], 9357408: 1, 9357440: 2, 9357472: 3, 28258000: [1, 2, 3, [1, 2, 3]]}
    #
    # key = 139907464678864
    # print(id(a) == key)  # True
    # print(id(b) == key)  # False
    # print(id(a) == id(memo[key]))  # False
    # print(id(b) == id(memo[key]))  # True
    #
    # in other words:
    # memo[id_of_initial_object] = copy_of_initial_object
    def __deepcopy__(self, memodict: dict = {}) -> object:
        import copy

        cls = self.__class__
        result = cls.__new__(cls)
        memodict[id(self)] = result
        for key, value in self.__dict__.items():
            setattr(result, key, copy.deepcopy(value, memodict))
        return result
