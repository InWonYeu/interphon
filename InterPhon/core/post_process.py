import numpy as np
from InterPhon.util import MatrixLike, AtomType, SelectIndex, FilePath, File, KptPath
from InterPhon.util import k_points, Symmetry2D
from InterPhon.core import UnitCell
from InterPhon.core import SuperCell
from InterPhon.core import PostArgument
from InterPhon.core import PreProcess
from InterPhon.inout import vasp, aims, espresso


class PostProcess(PreProcess):
    """
    Post process class to control post-process.
    This child class is inherited from the :class:`core.PreProcess` parent class.
    The information required in post-process is stored in the instance variables of this class.
    The instance variables are set by the :class:`core.PostProcess.set_user_arg`, :class:`core.PostProcess.set_reciprocal_lattice`,
    :class:`core.PostProcess.set_force_constant`, and :class:`core.PostProcess.set_k_points` methods.
    From the instance variables, eigen-frequency and eigen-mode are calculated using the :class:`core.PostProcess.eval_phonon` method,
    which will be further analyzed by employing the modules in analysis sub-package.

    :param in_file_unit_cell: Path of unit cell input file
    :type in_file_unit_cell: str
    :param in_file_super_cell: Path of super cell input file
    :type in_file_super_cell: str
    :param code_name: Specification of the file-format by a DFT program, defaults to vasp
    :type code_name: str
    :param user_arg: Instance of PostArgument class
    :type user_arg: :class:`core.PostArgument`
    :param unit_cell: Instance of UnitCell class
    :type unit_cell: :class:`core.UnitCell`
    :param super_cell: Instance of SuperCell class
    :type super_cell: :class:`core.SuperCell`
    """
    num_figure = 0

    def __init__(self, in_file_unit_cell: FilePath,
                 in_file_super_cell: FilePath,
                 code_name: str = 'vasp',
                 user_arg=PostArgument(),
                 unit_cell=UnitCell(),
                 super_cell=SuperCell()):
        """
        Constructor of PostProcess class.
        """
        super(PostProcess, self).__init__(user_arg,
                                          unit_cell,
                                          super_cell)

        self.unit_cell.read_unit_cell(in_file_unit_cell, code_name=code_name)
        self.unit_cell.set_mass_true()
        self.super_cell.read_unit_cell(in_file_super_cell, code_name=code_name)

        self.reciprocal_matrix = np.empty((3, 3))
        self.force_constant = np.empty((len(self.super_cell.atom_type) * 3,
                                        len(self.unit_cell.atom_true) * 3))
        self.k_points: KptPath = []
        self.auto_k_points = []
        self.dyn_matrix = np.empty((len(self.k_points),
                                    len(self.unit_cell.xyz_true),
                                    len(self.unit_cell.xyz_true)), dtype=complex)
        self.w_q = np.empty((len(self.k_points),
                             len(self.unit_cell.xyz_true)), dtype=float)
        self.v_q = np.empty((len(self.k_points),
                             len(self.unit_cell.xyz_true),
                             len(self.unit_cell.xyz_true)), dtype=complex)

    def set_user_arg(self, dict_args: dict) -> None:
        """
        Set a **PostArgument** instance from the information given by user.

        :param dict_args: Argument dictionary given by user
        :type dict_args: dict
        """
        super(PostProcess, self).set_user_arg(dict_args)
        self.user_arg.check_match_argument(self.unit_cell, self.super_cell)
        self.super_cell.set_super_ind_true(self.unit_cell, self.user_arg)

    def set_reciprocal_lattice(self) -> None:
        """
        Set the instance variable (**self.reciprocal_matrix**).
        """
        _volume = np.dot(self.unit_cell.lattice_matrix[0, 0:3],
                         np.cross(self.unit_cell.lattice_matrix[1, 0:3],
                                  self.unit_cell.lattice_matrix[2, 0:3]))
        for i in range(3):
            self.reciprocal_matrix[i, 0:3] = 2 * np.pi * np.cross(self.unit_cell.lattice_matrix[(i + 1) % 3, 0:3],
                                                                  self.unit_cell.lattice_matrix[(i + 2) % 3, 0:3]) / _volume

    def set_force_constant(self, force_files: FilePath,
                           code_name: str = 'vasp',
                           sym_flag: bool = True) -> None:
        """
        Set the instance variable (**self.force_constant**).

        :param force_files: Path of DFT output files which contain atomic forces
        :type force_files: str
        :param code_name: Specification of the file-format by a DFT program, defaults to vasp
        :type code_name: str
        :param sym_flag: Specify whether to use symmetry operation, defaults to `True`
        :type sym_flag: bool
        """
        _ind_pbc = self.user_arg.periodicity.nonzero()[0]
        if sym_flag:
            if _ind_pbc.shape[0] != 2:
                print('Caution:')
                print('Current version supports symmetry functionality only for 2D periodic systems.')
                print('"-sym" is changed from "{0}" to "False".'.format(sym_flag))
                sym_flag = False

        if code_name == 'vasp':
            if sym_flag:
                self.sym = Symmetry2D(self.unit_cell, self.super_cell, self.user_arg)
                _, _, _ = self.sym.search_point_group()
                _, _, _, _ = self.sym.search_image_atom()
                self.sym.search_self_image_atom()
                self.sym.search_independent_displacement()
                self.sym.gen_additional_displacement()

                force_ind = 0
                for i, require in enumerate(self.sym.require_atom):
                    num_of_calculation = 2 * (len(self.sym.independent_additional_displacement_cart[i]) + 1)

                    _independent_displace = [self.sym.independent_by_single_displacement_cart[i][_v] for _v in range(len(self.sym.independent_by_single_displacement_cart[i]))]
                    _additional_displace = self.sym.independent_additional_displacement_cart[i]
                    for __v in range(len(_additional_displace)):
                        _independent_displace.append(_additional_displace[__v])
                    to_cart_displace = np.transpose(np.array(_independent_displace))
                    to_random_displace = np.linalg.inv(to_cart_displace)

                    forward_force = np.empty((len(self.super_cell.atom_type) * 3, 3))
                    backward_force = np.empty((len(self.super_cell.atom_type) * 3, 3))
                    for j in range(num_of_calculation):
                        if force_ind % 2 == 0:
                            _forward_matrix = vasp.read_output_lines(force_files[force_ind], len(self.super_cell.atom_type))
                            __forward_matrix = _forward_matrix.copy()
                        elif force_ind % 2 == 1:
                            _backward_matrix = vasp.read_output_lines(force_files[force_ind], len(self.super_cell.atom_type))
                            __backward_matrix = _backward_matrix.copy()

                            if j // 2 == 0:
                                for k, (W_index, W_cart) in enumerate(zip(self.sym.independent_by_W_index[i], self.sym.independent_by_W_displacement_cart[i]), start=0):
                                    _forward_rot = W_cart @ np.transpose(_forward_matrix)
                                    _backward_rot = W_cart @ np.transpose(_backward_matrix)

                                    _image_sindex = [self.super_cell.atom_true[v] for i, v in enumerate(self.sym.same_supercell_index_select[W_index][require])]
                                    _original_sindex = [self.super_cell.atom_true[i] for i, v in enumerate( self.sym.same_supercell_index_select[W_index][require])]

                                    __forward_matrix[_image_sindex] = np.transpose(_forward_rot)[_original_sindex]
                                    __backward_matrix[_image_sindex] = np.transpose(_backward_rot)[_original_sindex]

                                    forward_force[:, j // 2 + k] = __forward_matrix.reshape([self.force_constant.shape[0], ])
                                    backward_force[:, j // 2 + k] = __backward_matrix.reshape([self.force_constant.shape[0], ])
                                    kk = k
                            else:
                                for _, _ in enumerate(self.sym.independent_additional_displacement_cart[i]):
                                    forward_force[:, j // 2 + kk] = __forward_matrix.reshape([self.force_constant.shape[0], ])
                                    backward_force[:, j // 2 + kk] = __backward_matrix.reshape([self.force_constant.shape[0], ])

                        force_ind += 1

                    _dif_force = - (forward_force - backward_force) @ to_random_displace / (2 * self.user_arg.displacement * 10 ** (-10))
                    self.force_constant[:, 3 * require: 3 * (require + 1)] = _dif_force

                _original_basis = np.transpose(self.unit_cell.lattice_matrix.copy())
                to_cart_coord = _original_basis / np.linalg.norm(_original_basis, axis=0)
                to_direct_coord = np.linalg.inv(to_cart_coord)
                for _point_group_ind, _not_require in zip(self.sym.point_group_ind, self.sym.not_require_atom):
                    W_in_cart = to_cart_coord @ self.sym.W_select[_point_group_ind] @ to_direct_coord

                    for _super_index, _super_same_index in enumerate(self.sym.same_supercell_index_select[_point_group_ind][_not_require]):
                        self.force_constant[3 * self.super_cell.atom_true[_super_index]: 3 * (self.super_cell.atom_true[_super_index] + 1),
                        3 * _not_require: 3 * (_not_require + 1)] \
                            = np.linalg.inv(W_in_cart) \
                              @ self.force_constant.copy()[3 * self.super_cell.atom_true[_super_same_index]: 3 * (self.super_cell.atom_true[_super_same_index] + 1),
                                3 * self.sym.same_index_select[_point_group_ind][0][_not_require]: 3 * (self.sym.same_index_select[_point_group_ind][0][_not_require] + 1)] \
                              @ W_in_cart

            else:
                for _ind_file, _force_file in enumerate(force_files):
                    if _ind_file % 2 == 0:
                        _forward_matrix = vasp.read_output_lines(_force_file, len(self.super_cell.atom_type))

                    elif _ind_file % 2 == 1:
                        _backward_matrix = vasp.read_output_lines(_force_file, len(self.super_cell.atom_type))

                        _dif_force = - (_forward_matrix - _backward_matrix) / (2 * self.user_arg.displacement * 10 ** (-10))
                        self.force_constant[:, _ind_file // 2] = _dif_force.reshape([self.force_constant.shape[0], ])

        elif code_name == 'espresso':
            if sym_flag:
                self.sym = Symmetry2D(self.unit_cell, self.super_cell, self.user_arg)
                _, _, _ = self.sym.search_point_group()
                _, _, _, _ = self.sym.search_image_atom()
                self.sym.search_self_image_atom()
                self.sym.search_independent_displacement()
                self.sym.gen_additional_displacement()

                force_ind = 0
                for i, require in enumerate(self.sym.require_atom):
                    num_of_calculation = 2 * (len(self.sym.independent_additional_displacement_cart[i]) + 1)

                    _independent_displace = [self.sym.independent_by_single_displacement_cart[i][_v] for _v in range(len(self.sym.independent_by_single_displacement_cart[i]))]
                    _additional_displace = self.sym.independent_additional_displacement_cart[i]
                    for __v in range(len(_additional_displace)):
                        _independent_displace.append(_additional_displace[__v])
                    to_cart_displace = np.transpose(np.array(_independent_displace))
                    to_random_displace = np.linalg.inv(to_cart_displace)

                    forward_force = np.empty((len(self.super_cell.atom_type) * 3, 3))
                    backward_force = np.empty((len(self.super_cell.atom_type) * 3, 3))
                    for j in range(num_of_calculation):
                        if force_ind % 2 == 0:
                            _forward_matrix = espresso.read_output_lines(force_files[force_ind], len(self.super_cell.atom_type))
                            __forward_matrix = _forward_matrix.copy()
                        elif force_ind % 2 == 1:
                            _backward_matrix = espresso.read_output_lines(force_files[force_ind], len(self.super_cell.atom_type))
                            __backward_matrix = _backward_matrix.copy()

                            if j // 2 == 0:
                                for k, (W_index, W_cart) in enumerate(zip(self.sym.independent_by_W_index[i], self.sym.independent_by_W_displacement_cart[i]), start=0):
                                    _forward_rot = W_cart @ np.transpose(_forward_matrix)
                                    _backward_rot = W_cart @ np.transpose(_backward_matrix)

                                    _image_sindex = [self.super_cell.atom_true[v] for i, v in enumerate(self.sym.same_supercell_index_select[W_index][require])]
                                    _original_sindex = [self.super_cell.atom_true[i] for i, v in enumerate(self.sym.same_supercell_index_select[W_index][require])]

                                    __forward_matrix[_image_sindex] = np.transpose(_forward_rot)[_original_sindex]
                                    __backward_matrix[_image_sindex] = np.transpose(_backward_rot)[_original_sindex]

                                    forward_force[:, j // 2 + k] = __forward_matrix.reshape([self.force_constant.shape[0], ])
                                    backward_force[:, j // 2 + k] = __backward_matrix.reshape([self.force_constant.shape[0], ])
                                    kk = k
                            else:
                                for _, _ in enumerate(self.sym.independent_additional_displacement_cart[i]):
                                    forward_force[:, j // 2 + kk] = __forward_matrix.reshape([self.force_constant.shape[0], ])
                                    backward_force[:, j // 2 + kk] = __backward_matrix.reshape([self.force_constant.shape[0], ])

                        force_ind += 1

                    _dif_force = - (forward_force - backward_force) @ to_random_displace / (2 * self.user_arg.displacement * 10 ** (-10))
                    self.force_constant[:, 3 * require: 3 * (require + 1)] = _dif_force

                _original_basis = np.transpose(self.unit_cell.lattice_matrix.copy())
                to_cart_coord = _original_basis / np.linalg.norm(_original_basis, axis=0)
                to_direct_coord = np.linalg.inv(to_cart_coord)
                for _point_group_ind, _not_require in zip(self.sym.point_group_ind, self.sym.not_require_atom):
                    W_in_cart = to_cart_coord @ self.sym.W_select[_point_group_ind] @ to_direct_coord

                    for _super_index, _super_same_index in enumerate(self.sym.same_supercell_index_select[_point_group_ind][_not_require]):
                        self.force_constant[3 * self.super_cell.atom_true[_super_index]: 3 * (self.super_cell.atom_true[_super_index] + 1),
                        3 * _not_require: 3 * (_not_require + 1)] \
                            = np.linalg.inv(W_in_cart) \
                              @ self.force_constant.copy()[3 * self.super_cell.atom_true[_super_same_index]: 3 * (self.super_cell.atom_true[_super_same_index] + 1),
                                3 * self.sym.same_index_select[_point_group_ind][0][_not_require]: 3 * (self.sym.same_index_select[_point_group_ind][0][_not_require] + 1)] \
                              @ W_in_cart

            else:
                for _ind_file, _force_file in enumerate(force_files):
                    if _ind_file % 2 == 0:
                        _forward_matrix = espresso.read_output_lines(_force_file, len(self.super_cell.atom_type))

                    elif _ind_file % 2 == 1:
                        _backward_matrix = espresso.read_output_lines(_force_file, len(self.super_cell.atom_type))

                        _dif_force = - (_forward_matrix - _backward_matrix) / (2 * self.user_arg.displacement * 10 ** (-10))
                        self.force_constant[:, _ind_file // 2] = _dif_force.reshape([self.force_constant.shape[0], ])

        elif code_name == 'aims':
            if sym_flag:
                self.sym = Symmetry2D(self.unit_cell, self.super_cell, self.user_arg)
                _, _, _ = self.sym.search_point_group()
                _, _, _, _ = self.sym.search_image_atom()
                self.sym.search_self_image_atom()
                self.sym.search_independent_displacement()
                self.sym.gen_additional_displacement()

                force_ind = 0
                for i, require in enumerate(self.sym.require_atom):
                    num_of_calculation = 2 * (len(self.sym.independent_additional_displacement_cart[i]) + 1)

                    _independent_displace = [self.sym.independent_by_single_displacement_cart[i][_v] for _v in range(len(self.sym.independent_by_single_displacement_cart[i]))]
                    _additional_displace = self.sym.independent_additional_displacement_cart[i]
                    for __v in range(len(_additional_displace)):
                        _independent_displace.append(_additional_displace[__v])
                    to_cart_displace = np.transpose(np.array(_independent_displace))
                    to_random_displace = np.linalg.inv(to_cart_displace)

                    forward_force = np.empty((len(self.super_cell.atom_type) * 3, 3))
                    backward_force = np.empty((len(self.super_cell.atom_type) * 3, 3))
                    for j in range(num_of_calculation):
                        if force_ind % 2 == 0:
                            _forward_matrix = aims.read_output_lines(force_files[force_ind], len(self.super_cell.atom_type))
                            __forward_matrix = _forward_matrix.copy()
                        elif force_ind % 2 == 1:
                            _backward_matrix = aims.read_output_lines(force_files[force_ind], len(self.super_cell.atom_type))
                            __backward_matrix = _backward_matrix.copy()

                            if j // 2 == 0:
                                for k, (W_index, W_cart) in enumerate(zip(self.sym.independent_by_W_index[i], self.sym.independent_by_W_displacement_cart[i]), start=0):
                                    _forward_rot = W_cart @ np.transpose(_forward_matrix)
                                    _backward_rot = W_cart @ np.transpose(_backward_matrix)

                                    _image_sindex = [self.super_cell.atom_true[v] for i, v in enumerate(self.sym.same_supercell_index_select[W_index][require])]
                                    _original_sindex = [self.super_cell.atom_true[i] for i, v in enumerate(self.sym.same_supercell_index_select[W_index][require])]

                                    __forward_matrix[_image_sindex] = np.transpose(_forward_rot)[_original_sindex]
                                    __backward_matrix[_image_sindex] = np.transpose(_backward_rot)[_original_sindex]

                                    forward_force[:, j // 2 + k] = __forward_matrix.reshape([self.force_constant.shape[0], ])
                                    backward_force[:, j // 2 + k] = __backward_matrix.reshape([self.force_constant.shape[0], ])
                                    kk = k
                            else:
                                for _, _ in enumerate(self.sym.independent_additional_displacement_cart[i]):
                                    forward_force[:, j // 2 + kk] = __forward_matrix.reshape([self.force_constant.shape[0], ])
                                    backward_force[:, j // 2 + kk] = __backward_matrix.reshape([self.force_constant.shape[0], ])

                        force_ind += 1

                    _dif_force = - (forward_force - backward_force) @ to_random_displace / (2 * self.user_arg.displacement * 10 ** (-10))
                    self.force_constant[:, 3 * require: 3 * (require + 1)] = _dif_force

                _original_basis = np.transpose(self.unit_cell.lattice_matrix.copy())
                to_cart_coord = _original_basis / np.linalg.norm(_original_basis, axis=0)
                to_direct_coord = np.linalg.inv(to_cart_coord)
                for _point_group_ind, _not_require in zip(self.sym.point_group_ind, self.sym.not_require_atom):
                    W_in_cart = to_cart_coord @ self.sym.W_select[_point_group_ind] @ to_direct_coord

                    for _super_index, _super_same_index in enumerate(self.sym.same_supercell_index_select[_point_group_ind][_not_require]):
                        self.force_constant[3 * self.super_cell.atom_true[_super_index]: 3 * (self.super_cell.atom_true[_super_index] + 1),
                        3 * _not_require: 3 * (_not_require + 1)] \
                            = np.linalg.inv(W_in_cart) \
                              @ self.force_constant.copy()[3 * self.super_cell.atom_true[_super_same_index]: 3 * (self.super_cell.atom_true[_super_same_index] + 1),
                                3 * self.sym.same_index_select[_point_group_ind][0][_not_require]: 3 * (self.sym.same_index_select[_point_group_ind][0][_not_require] + 1)] \
                              @ W_in_cart

            else:
                for _ind_file, _force_file in enumerate(force_files):
                    if _ind_file % 2 == 0:
                        _forward_matrix = aims.read_output_lines(_force_file, len(self.super_cell.atom_type))

                    elif _ind_file % 2 == 1:
                        _backward_matrix = aims.read_output_lines(_force_file, len(self.super_cell.atom_type))

                        _dif_force = - (_forward_matrix - _backward_matrix) / (2 * self.user_arg.displacement * 10 ** (-10))
                        self.force_constant[:, _ind_file // 2] = _dif_force.reshape([self.force_constant.shape[0], ])

    def set_k_points(self, k_file: FilePath) -> None:
        """
        Set the instance variable (**self.k_points**) by reading **KPOINTS** file given in VASP format.
        The following four samplings of the wave vector k are available:
        **1) Gamma-centered**
        **2) Monkhorst-Pack**
        **3) Line-path generation** (for band plot)
        **4) Explicit designation**

        For more details, see the following reference:
        "The Basics of Electronic Structure Theory for Periodic Systems, Frontiers in chemistry 7, 1 (2019)".

        :param k_file: Path of KPOINTS file
        :type k_file: str
        """
        try:
            with open(k_file, 'r') as infile:
                lines = infile.readlines()
        except IOError:
            print("\nFail to open '{0}' file".format(k_file))
            print("Please check the path of k-points file\n")
            raise

        self.k_points: KptPath = []
        self.auto_k_points = []

        _ind_pbc = self.user_arg.periodicity.nonzero()[0]

        if lines[2].split()[0][0] in ('G', 'g'):
            # Gamma-centered grid of k-points
            self.k_points, self.auto_k_points = k_points.gamma_centered(lines, _ind_pbc)

        elif lines[2].split()[0][0] in ('M', 'm'):
            # Monkhorst-Pack grid of k-points
            self.k_points, self.auto_k_points = k_points.monkhorst_pack(lines, _ind_pbc)

        elif lines[2].split()[0][0] in ('L', 'l'):
            # Line-path of k-points
            self.k_points = k_points.line_path(lines, _ind_pbc)

        elif lines[2].split()[0][0] in ('R', 'r'):
            # Explicit list of k-points
            self.k_points = k_points.explicit_reciprocal(lines)

    def eval_phonon(self) -> None:
        """
        Set the instance variables, eigen-frequency (**self.w_q**) and corresponding eigen-mode (**self.v_q**).
        """
        if len(self.k_points) != 0:
            self.dyn_matrix = np.empty((len(self.k_points),
                                        len(self.unit_cell.xyz_true),
                                        len(self.unit_cell.xyz_true)), dtype=complex)
            self.w_q = np.empty((len(self.k_points),
                                 len(self.unit_cell.xyz_true)), dtype=float)
            self.v_q = np.empty((len(self.k_points),
                                 len(self.unit_cell.xyz_true),
                                 len(self.unit_cell.xyz_true)), dtype=complex)

        _ind_k = 0
        _enlarge = 1
        for ind, value in enumerate(self.user_arg.periodicity):
            if value:
                _enlarge = _enlarge * self.user_arg.enlargement[ind]
        _ind_pbc = self.user_arg.periodicity.nonzero()[0]

        for k_point in self.k_points:
            q = np.dot(k_point, self.reciprocal_matrix)
            _dyn_matrix = np.zeros([len(self.unit_cell.xyz_true), len(self.unit_cell.xyz_true)], dtype=complex)

            for s1, satom_ind in enumerate(self.super_cell.xyz_true):
                for s2, atom_ind in enumerate(self.unit_cell.xyz_true):
                    pos_vector = self.super_cell.atom_cart[satom_ind, 0:3] - self.unit_cell.atom_cart[atom_ind, 0:3]

                    if _ind_pbc.shape[0] == 0:
                        pass

                    elif _ind_pbc.shape[0] == 1:
                        for first in (-1, 0, 1):
                            _tmp_pos = pos_vector + first * self.super_cell.lattice_matrix[_ind_pbc[0], 0:3]
                            if np.dot(pos_vector, pos_vector) > np.dot(_tmp_pos, _tmp_pos):
                                pos_vector = _tmp_pos.copy()

                    elif _ind_pbc.shape[0] == 2:
                        for first in (-1, 0, 1):
                            for second in (-1, 0, 1):
                                _tmp_pos = pos_vector + first * self.super_cell.lattice_matrix[_ind_pbc[0], 0:3] \
                                           + second * self.super_cell.lattice_matrix[_ind_pbc[1], 0:3]
                                if np.dot(pos_vector, pos_vector) > np.dot(_tmp_pos, _tmp_pos):
                                    pos_vector = _tmp_pos.copy()

                    elif _ind_pbc.shape[0] == 3:
                        for first in (-1, 0, 1):
                            for second in (-1, 0, 1):
                                for third in (-1, 0, 1):
                                    _tmp_pos = pos_vector + first * self.super_cell.lattice_matrix[_ind_pbc[0], 0:3] \
                                               + second * self.super_cell.lattice_matrix[_ind_pbc[1], 0:3] \
                                               + third * self.super_cell.lattice_matrix[_ind_pbc[2], 0:3]
                                    if np.dot(pos_vector, pos_vector) > np.dot(_tmp_pos, _tmp_pos):
                                        pos_vector = _tmp_pos.copy()

                    _dyn_matrix[(s1 // (_enlarge * 3)) * 3 + (s1 % 3), s2] = \
                        _dyn_matrix[(s1 // (_enlarge * 3)) * 3 + (s1 % 3), s2] \
                        + complex(self.force_constant[satom_ind * 3 + (s1 % 3), s2], 0) \
                        * np.exp(1j * complex(np.dot(q, pos_vector), 0))

            _mass_true = np.sqrt(np.dot(self.unit_cell.mass_true.reshape([len(self.unit_cell.xyz_true), 1]),
                                        self.unit_cell.mass_true.reshape([1, len(self.unit_cell.xyz_true)])))
            self.dyn_matrix[_ind_k, :, :] = _dyn_matrix / _mass_true
            _eig_w, self.v_q[_ind_k, :, :] = np.linalg.eig(self.dyn_matrix[_ind_k, :, :])

            self.w_q[_ind_k, :] = (np.sqrt(_eig_w).real - np.abs(np.sqrt(_eig_w).imag)) / (2 * np.pi) / 10 ** 12  # THz
            self.w_q[_ind_k, :] = self.w_q[_ind_k, np.argsort(self.w_q[_ind_k, :])]
            self.v_q[_ind_k, :, :] = self.v_q[_ind_k, np.argsort(self.w_q[_ind_k, :]), :]
            _ind_k = _ind_k + 1

        if _ind_pbc.shape[0] == 0:
            with open('./freqency_at_gamma_point.dat', 'w') as outfile:
                comment = "Discrete frequency of non-periodic system"
                outfile.write("%s" % comment + '\n')
                outfile.write("%s" % '    Frequency (THz)' + '\n')

                for _freq in self.w_q[0]:
                    line = ' %16.9f ' % _freq
                    outfile.write("%s" % line + '\n')

    # def write_dos(self, out_file: FilePath = 'total_dos.dat', sigma: float = 0.0, num_dos: int = 200,
    #               partial_dos: bool = False, plot: bool = False) -> File:
    #     """
    #     Method of PostProcess class.
    #     Process to write and plot phonon DOS.
    #
    #     usage:
    #     " >>> instance_of_PostProcess.write_dos(out_file='total_dos.dat', sigma=0.1, num_dos=1000,
    #     partial_dos=True, plot=True)"
    #
    #     :param out_file: (str) Name of output file for total dos.
    #     :param sigma: (float) Value of sigma for gaussian smearing,
    #     if the sigma value is 0.0, DOS will be presented by the linear tetrahedron method.
    #     :param num_dos: (int) The number of energy points evaluated for DOS.
    #     :param partial_dos: (bool) Write (True) output file for partial dos or not (False).
    #     :param plot: (bool) Plot (True) phonon DOS or not (False).
    #     :return: (File)
    #     """
    #     minimum_freq = np.min(self.w_q)
    #     maximum_freq = np.max(self.w_q)
    #     _freq = np.arange(minimum_freq - 2, maximum_freq + 2, (maximum_freq - minimum_freq + 4) / num_dos)
    #     _pdos = np.zeros((len(self.unit_cell.xyz_true), _freq.shape[0]))
    #
    #     if sigma == 0.0:
    #         # Linear Tetrahedron Method for Brillouin zone integration
    #         _ind_pbc = self.user_arg.periodicity.nonzero()[0]
    #         if _ind_pbc.shape[0] == 0:
    #             _freq = self.w_q[0, :]
    #             _tdos = np.ones(_freq.shape[0])
    #
    #         elif _ind_pbc.shape[0] == 1:
    #             _pdos = tetrahedron_1d(_freq, _pdos, self.k_points, self.w_q, self.v_q, self.auto_k_points, _ind_pbc)
    #             _tdos = _pdos.sum(axis=0)
    #
    #         elif _ind_pbc.shape[0] == 2:
    #             _pdos = tetrahedron_2d(_freq, _pdos, self.k_points, self.w_q, self.v_q, self.auto_k_points, _ind_pbc)
    #             _tdos = _pdos.sum(axis=0)
    #
    #         elif _ind_pbc.shape[0] == 3:
    #             _pdos = tetrahedron_3d(_freq, _pdos, self.k_points, self.w_q, self.v_q, self.auto_k_points, _ind_pbc)
    #             _tdos = _pdos.sum(axis=0)
    #
    #     else:
    #         # Gaussian Smearing Method for Brillouin zone integration
    #         for eig_freqs, eig_modes in zip(self.w_q, self.v_q):
    #             for ind_freq, eig_freq in enumerate(eig_freqs):
    #                 for ind_mode, _ in enumerate(self.unit_cell.xyz_true):
    #                     _pdos[ind_mode, :] = _pdos[ind_mode, :] \
    #                                     + 1 / (sigma * np.sqrt(2 * np.pi)) \
    #                                     * np.exp(- (_freq - eig_freq) ** 2 / (2 * sigma ** 2)) \
    #                                     / len(self.k_points) \
    #                                     * (abs(eig_modes[ind_freq, ind_mode]) ** 2)
    #
    #         _tdos = _pdos.sum(axis=0)
    #
    #     with open(out_file + '/total_dos.dat', 'w') as outfile:
    #         if sigma == 0.0:
    #             comment = "Total phonon DOS by Linear Tetrahedron Method"
    #             outfile.write("%s" % comment + '\n')
    #
    #         elif sigma != 0.0:
    #             comment = "Total phonon DOS by Gaussian Smearing with Sigma = %f" % sigma
    #             outfile.write("%s" % comment + '\n')
    #
    #         outfile.write("%s" %
    #                       '    Frequency (THz)' +
    #                       '    Total_Density_of_State' + '\n')
    #
    #         for x, y in zip(_freq, _tdos):
    #             line = ' %16.9f ' % x + ' %16.9f ' % y
    #             outfile.write("%s" % line + '\n')
    #
    #     if partial_dos is True:
    #         with open(out_file + '/partial_dos.dat', 'w') as outfile:
    #             if sigma == 0.0:
    #                 comment = "Partial phonon DOS by Linear Tetrahedron Method"
    #                 outfile.write("%s" % comment + '\n')
    #
    #             elif sigma != 0.0:
    #                 comment = "Partial phonon DOS by Gaussian Smearing with Sigma = %f" % sigma
    #                 outfile.write("%s" % comment + '\n')
    #
    #             line = "%s" % '    Frequency (THz)'
    #             for xyz_true in self.unit_cell.xyz_true:
    #                 line += '   pdos-atom-{0:0>3}  '.format(xyz_true)
    #             outfile.write(line + '\n')
    #
    #             for x_ind, x in enumerate(_freq):
    #                 line = ' %16.9f ' % x
    #                 for ind, _ in enumerate(self.unit_cell.xyz_true):
    #                     line += ' %16.9f ' % _pdos[ind, x_ind]
    #                 outfile.write(line + '\n')
    #
    #         if plot is True:
    #             from matplotlib import pyplot as plt
    #             self.num_figure = self.num_figure + 1
    #             plt.figure(self.num_figure)
    #             plt.plot(_freq, _tdos, color='black', label='total dos')
    #             for ind, xyz_true in enumerate(self.unit_cell.xyz_true):
    #                 plt.plot(_freq, _pdos[ind], linestyle='--', label='atom-{0:0>4}'.format(xyz_true))
    #             plt.xlabel('Frequency (THz)')
    #             plt.ylabel('DOS (a.u.)')
    #             plt.legend(fontsize='x-large')
    #             plt.show()
    #
    #     elif plot is True:
    #         from matplotlib import pyplot as plt
    #         self.num_figure = self.num_figure + 1
    #         plt.figure(self.num_figure)
    #         plt.plot(_freq, _tdos, color='black', label='total dos')
    #         plt.xlabel('Frequency (THz)')
    #         plt.ylabel('DOS (a.u.)')
    #         plt.legend(fontsize='x-large')
    #         plt.show()
    #
    # def write_band(self, out_file: FilePath = 'band.dat', plot: bool = False) -> File:
    #     """
    #     Method of PostProcess class.
    #     Process to write and plot phonon band.
    #
    #     usage:
    #     " >>> instance_of_PostProcess.write_band(out_file='band.dat', plot=True)"
    #
    #     :param out_file: (str) Name of output file for phonon band.
    #     :param plot: (bool) Plot (True) phonon band or not (False).
    #     :return: (File)
    #     """
    #
    #     with open(out_file, 'w') as outfile:
    #         comment = "Phonon Band"
    #         outfile.write("%s" % comment + '\n')
    #         outfile.write("%s" %
    #                       '    K_Points_Path' +
    #                       '    Frequency (THz)' + '\n')
    #
    #         k_points_length = np.zeros(len(self.k_points))
    #         for ind, kpt in enumerate(self.k_points[1:], 1):
    #             k_points_length[ind] = k_points_length[ind - 1] \
    #                                    + np.sqrt(np.dot(self.k_points[ind] - self.k_points[ind - 1],
    #                                                     self.k_points[ind] - self.k_points[ind - 1]))
    #
    #         for x, y_set in zip(k_points_length, self.w_q):
    #             line = ' %16.9f ' % x
    #             for y in y_set:
    #                 line = line + ' %16.9f ' % y
    #             outfile.write("%s" % line + '\n')
    #
    #     if plot is True:
    #         from matplotlib import pyplot as plt
    #         self.num_figure = self.num_figure + 1
    #         plt.figure(self.num_figure)
    #         plt.plot(k_points_length, self.w_q, color='black')
    #         plt.xlabel('K points')
    #         plt.ylabel('Frequency (THz)')
    #         plt.show()
    #
    # def write_mode(self, out_file: FilePath = 'mode',
    #                mode_inds: SelectIndex = [0],
    #                k_point: MatrixLike = np.array([0.0, 0.0, 0.0]),
    #                num_images: int = 30,
    #                plot: bool = False) -> File:
    #     """
    #     Method of PostProcess class.
    #     Process to write and plot phonon mode.
    #     Plotting phonon mode is supported with the help of Atomic Simulation Environment (ASE) package.
    #     (see, https://wiki.fysik.dtu.dk/ase/index.html)
    #     Therefore, in order to see phonon mode, the ASE package must be installed and
    #     the position of the package has to be automatically found by the Python interpreter (plot=True).
    #
    #     usage:
    #     " >>> instance_of_PostProcess.write_mode(out_file='mode', mode_inds=[0,1],
    #     k_point=[0.0, 0.0, 0.0], num_images=10, plot=True)"
    #
    #     :param out_file: (str) Name of output file for phonon mode.
    #     :param mode_inds: (List[int]) Mode index where '0' indicates the mode of the lowest frequency.
    #     :param k_point: (np.ndarray[float]) '(3,) size' k-point where corresponding phonon mode will be evaluated.
    #     :param num_images: (int) The number of images to see the oscillating phonon mode.
    #     :param plot: (bool) Plot (True) phonon mode or not (False).
    #     :return: (File)
    #     """
    #     k_point = np.array(k_point)
    #
    #     _ind_k = None
    #     for ind, kpt in enumerate(self.k_points):
    #         if np.allclose(k_point, kpt):
    #             _ind_k = ind
    #     if _ind_k is None:
    #         raise error.Not_Specified_Kpath_Error
    #
    #     mode = self.v_q[_ind_k, :, :]
    #
    #     _current_position_true = self.unit_cell.atom_cart.copy()[self.unit_cell.atom_true]
    #     for mode_ind in mode_inds:
    #         for ind, x in enumerate(np.linspace(0, 2 * np.pi, num_images, endpoint=False), 1):
    #             self.unit_cell.atom_cart[self.unit_cell.atom_true] = _current_position_true + \
    #                                                             np.sin(x) * mode[mode_ind, :].reshape((-1, 3)).real
    #
    #             if ind == 1:
    #                 lines = vasp.write_input_lines(self.unit_cell, "unknown system")
    #                 lines[7] = "Cartesian configuration= %4d" % ind + '\n'
    #             else:
    #                 line = vasp.write_input_lines(self.unit_cell, "unknown system")
    #                 line[7] = "Cartesian configuration= %4d" % ind + '\n'
    #                 lines.extend(line[7:])
    #         with open("XDATCAR" + "_%s_%d" % (out_file, mode_ind), 'w') as outfile:
    #             outfile.write("%s" % "".join(lines))
    #
    #     if plot is True:
    #         try:
    #             from ase.io.trajectory import Trajectory
    #             from ase.io import read, iread
    #             from ase.visualize import view
    #         except ImportError:
    #             print("\nThe parent directory of ase package must be included in 'sys.path'")
    #             print("Example: 'sys.path.append('C:\\Users\\user\\Documents')")
    #             print("if the parent directory of ase package is 'C:\\Users\\user\\Documents'\n")
    #             raise
    #
    #         atom = read('POSCAR')
    #         _current_position_true = atom.positions.copy()[self.unit_cell.atom_true]
    #         for mode_ind in mode_inds:
    #             traj = Trajectory("%s.%d.traj" % (out_file, mode_ind), 'w')
    #             for ind, x in enumerate(np.linspace(0, 2 * np.pi, num_images, endpoint=False), 1):
    #                 atom.positions[self.unit_cell.atom_true] = _current_position_true + \
    #                                                            np.sin(x) * mode[mode_ind, :].reshape((-1, 3)).real
    #                 traj.write(atom)
    #             traj.close()
    #             atoms = iread("%s.%d.traj" % (out_file, mode_ind))
    #             view(atoms)
    #
    # def write_thermal_properties(self, out_file: FilePath = 'thermal_properties.dat',
    #                              temp: MatrixLike = np.arange(0, 1001, 10),
    #                              plot: bool = False) -> File:
    #     """
    #     Method of PostProcess class.
    #     Process to write and plot thermal properties.
    #
    #     usage:
    #     " >>> instance_of_PostProcess.write_thermal_properties(out_file='thermal_properties.dat',
    #     temp=np.arange(0, 1201, 10), plot=True)"
    #
    #     :param out_file: (str) Name of output file for thermal properties.
    #     :param temp: (np.ndarray) Temperature array where thermal properties will be evaluated.
    #     :param plot: (bool) Plot (True) thermal properties or not (False).
    #     :return: (File)
    #     """
    #     with open(out_file, 'w') as outfile:
    #         comment = "Thermal Properties / atom"
    #         outfile.write("%s" % comment + '\n')
    #         outfile.write("%s" %
    #                       '    Temperature (K)' +
    #                       '    Free_energy (eV/atom)' +
    #                       '    Entropy (meV/K/atom)' +
    #                       '    Heat_capacity (meV/K/atom)' + '\n')
    #
    #         kb = 1.38 * 10 ** (-23) / (1.602 * 10 ** (-19))
    #         h = 6.626 * 10 ** (-34) / (1.602 * 10 ** (-19))
    #
    #         temp = np.array(temp)
    #         free_energy = np.zeros(temp.shape)
    #         entropy = np.zeros(temp.shape)
    #         heat_capacity = np.zeros(temp.shape)
    #
    #         for eig_freqs in self.w_q:
    #             for eig_freq in eig_freqs:
    #                 for ind, t in enumerate(temp):
    #                     if t < 10**(-6):
    #                         free_energy[ind] = free_energy[ind] + (1.0 / 2.0 * h * eig_freq) / len(self.k_points) / (
    #                                     self.w_q.shape[1] / 3)
    #                     else:
    #                         free_energy[ind] = free_energy[ind] + (1.0 / 2.0 * h * eig_freq
    #                                              + kb * t * np.log(1 - np.exp(-h * eig_freq / (kb * t)))) / \
    #                               len(self.k_points) / (self.w_q.shape[1] / 3)
    #                         entropy[ind] = entropy[ind] + (-kb * np.log(1 - np.exp(-h * eig_freq / (kb * t)))
    #                                      + 1 / t * h * eig_freq / (np.exp(h * eig_freq / (kb * t)) - 1)) / \
    #                           len(self.k_points) / (self.w_q.shape[1] / 3)
    #                         heat_capacity[ind] = heat_capacity[ind] + kb * np.power(h * eig_freq / (kb * t), 2) \
    #                                 * np.exp(h * eig_freq / (kb * t)) / \
    #                                 np.power(np.exp(h * eig_freq / (kb * t)) - 1, 2) / \
    #                                 len(self.k_points) / (self.w_q.shape[1] / 3)
    #
    #         for x, y1, y2, y3 in zip(temp, free_energy, entropy, heat_capacity):
    #             line = ' %20.9f ' % x + ' %20.9f ' % y1 + ' %20.9f ' % (y2 * 1000) + ' %20.9f ' % (y3 * 1000)
    #             outfile.write("%s" % line + '\n')
    #
    #     if plot is True:
    #         from matplotlib import pyplot as plt
    #         self.num_figure = self.num_figure + 1
    #         plt.figure(self.num_figure)
    #         plt.plot(temp, free_energy, label='Free energy (eV/atom)')
    #         plt.plot(temp, entropy * 1000, label='Entropy (meV/K/atom)')
    #         plt.plot(temp, heat_capacity * 1000, label='Heat capacity (meV/K/atom)')
    #         plt.xlabel('Temperature (K)')
    #         plt.ylabel('Thermal Properties')
    #         plt.legend(fontsize='x-large')
    #         plt.show()
