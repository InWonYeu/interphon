import numpy as np

# Searching lattice point group operations
W_candidate = [np.array([[0, 1],
                         [1, 0]]), np.array([[0, 1],
                                             [-1, 0]]), np.array([[0, -1],
                                                                  [1, 0]]), np.array([[0, -1],
                                                                                      [-1, 0]]),
               np.array([[0, 1],
                         [1, 1]]), np.array([[0, 1],
                                             [1, -1]]), np.array([[0, 1],
                                                                  [-1, 1]]), np.array([[0, -1],
                                                                                       [1, 1]]), np.array([[0, -1],
                                                                                                           [-1, 1]]), np.array([[0, -1],
                                                                                                                                [1, -1]]), np.array([[0, 1],
                                                                                                                                                     [-1, -1]]), np.array([[0, -1],
                                                                                                                                                                           [-1, -1]]),
               np.array([[1, 0],
                         [0, 1]]), np.array([[1, 0],
                                             [0, -1]]), np.array([[-1, 0],
                                                                  [0, 1]]), np.array([[-1, 0],
                                                                                      [0, -1]]),
               np.array([[1, 0],
                         [1, 1]]), np.array([[1, 0],
                                             [1, -1]]), np.array([[1, 0],
                                                                  [-1, 1]]), np.array([[-1, 0],
                                                                                       [1, 1]]), np.array([[-1, 0],
                                                                                                           [-1, 1]]), np.array([[-1, 0],
                                                                                                                                [1, -1]]), np.array([[1, 0],
                                                                                                                                                     [-1, -1]]), np.array([[-1, 0],
                                                                                                                                                                           [-1, -1]]),
               np.array([[1, 1],
                         [0, 1]]), np.array([[1, 1],
                                             [0, -1]]), np.array([[1, -1],
                                                                  [0, 1]]), np.array([[-1, 1],
                                                                                      [0, 1]]), np.array([[-1, -1],
                                                                                                          [0, 1]]), np.array([[-1, 1],
                                                                                                                              [0, -1]]), np.array([[1, -1],
                                                                                                                                                   [0, -1]]), np.array([[-1, -1],
                                                                                                                                                                        [0, -1]]),
               np.array([[1, 1],
                         [1, 0]]), np.array([[1, 1],
                                             [-1, 0]]), np.array([[1, -1],
                                                                  [1, 0]]), np.array([[-1, 1],
                                                                                      [1, 0]]), np.array([[-1, -1],
                                                                                                          [1, 0]]), np.array([[-1, 1],
                                                                                                                              [-1, 0]]), np.array([[1, -1],
                                                                                                                                                   [-1, 0]]), np.array([[-1, -1],
                                                                                                                                                                        [-1, 0]]),
               ]


class Symmetry2D(object):
    def __init__(self, unit_cell, super_cell, user_arg):
        self.unit_cell = unit_cell
        self.super_cell = super_cell
        self.user_arg = user_arg

        self.W_select = []
        self.w_select = []
        self.same_index_select = []
        self.point_group = None

        self.point_group_ind = []
        self.require_atom = []
        self.not_require_atom = []
        self.same_supercell_index_select = []

        self.point_group_for_self_require_atom = []

        self.independent_by_W_index = []
        self.independent_by_W_displacement_cart = []
        self.independent_by_single_displacement_cart = []

        self.independent_additional_displacement_cart = []

    # @property
    # def user_arg(self):
    #     return self.__user_arg
    #
    # @user_arg.setter
    # def user_arg(self, _user_arg):
    #     if isinstance(_user_arg, PreArgument):
    #         self.__user_arg = _user_arg
    #     else:
    #         ValueError(
    #             "'{0}' should be the instance of <class 'InterPhon.core.pre_check.PreArgument'>".format(_user_arg))
    #
    # @property
    # def unit_cell(self):
    #     return self.__unit_cell
    #
    # @unit_cell.setter
    # def unit_cell(self, _unit_cell):
    #     if isinstance(_unit_cell, UnitCell):
    #         self.__unit_cell = _unit_cell
    #     else:
    #         ValueError(
    #             "'{0}' should be the instance of <class 'InterPhon.core.unit_cell.UnitCell'>".format(_unit_cell))
    #
    # @property
    # def super_cell(self):
    #     return self.__super_cell
    #
    # @super_cell.setter
    # def super_cell(self, _super_cell):
    #     if isinstance(_super_cell, SuperCell):
    #         self.__super_cell = _super_cell
    #     else:
    #         ValueError(
    #             "'{0}' should be the instance of <class 'InterPhon.core.super_cell.SuperCell'>".format(_super_cell))

    def search_point_group(self):
        # metric tensor
        G_metric = np.dot(self.unit_cell.lattice_matrix.copy()[np.ix_(self.user_arg.periodicity.nonzero()[0],
                                                                      self.user_arg.periodicity.nonzero()[0])],
                          np.transpose(self.unit_cell.lattice_matrix.copy()[np.ix_(self.user_arg.periodicity.nonzero()[0],
                                                                                   self.user_arg.periodicity.nonzero()[0])]))

        rot_ind = []
        for ind, rot in enumerate(W_candidate):
            G_rotate = np.dot(np.transpose(rot), np.dot(G_metric, rot))
            if np.allclose(G_metric, G_rotate, atol=1e-06):
                rot_ind.append(ind)

        atom_true_original = np.transpose(self.unit_cell.atom_direct.copy()[self.unit_cell.atom_true, :])
        w_for_given_rot = []
        same_index = []
        for ind in rot_ind:
            _W_candidate = np.identity(3)
            _W_candidate[np.ix_(self.user_arg.periodicity.nonzero()[0],
                                self.user_arg.periodicity.nonzero()[0])] = W_candidate.copy()[ind]
            atom_true_rot = np.dot(_W_candidate, atom_true_original)

            # space group search
            w_candidate = [np.array([0.0, 0.0, 0.0])]
            # for index, value in enumerate(unit_cell.atom_true):
            #     other_atom_ind = [i for i in unit_cell.atom_true]
            #     other_atom_ind.remove(value)
            #     for other_index, other_ind in enumerate(other_atom_ind):
            #         if unit_cell.atom_type[other_ind] == unit_cell.atom_type[value]:
            #             w = np.round_(atom_true_original[:, other_index] - atom_true_rot[:, index], 6)
            #             if w[2] == 0.0:
            #                 w_, _ = np.modf(w)
            #                 if not np.allclose(w_, np.zeros([3, ]), atol=1e-06):
            #                     w_candidate.append(w_)
            # print('w_candidate=', w_candidate)
            # print(np.array(w_candidate).shape)

            trans_for_given_rot = []
            _same_index = []
            for w in w_candidate:
                atom_transform = atom_true_rot + w.reshape([3, 1])

                __same_index = []
                for index, value in enumerate(self.unit_cell.atom_true):
                    same_atom_type = [ind_ for ind_, val_ in enumerate(self.unit_cell.atom_true)
                                      if self.unit_cell.atom_type[val_] == self.unit_cell.atom_type[value]]

                    for _, same_atom_index in enumerate(same_atom_type):
                        delta_x = atom_transform[:, index] - atom_true_original[:, same_atom_index]  # atom-to-atom compare
                        delta_x_cart = np.matmul(np.transpose(self.unit_cell.lattice_matrix.copy()), delta_x - np.rint(delta_x))

                        if np.allclose(delta_x_cart, np.zeros([3, ]), atol=1e-06):
                            # if same_atom_index not in __same_index:
                            __same_index.append(same_atom_index)
                            # break

                    if len(__same_index) == len(self.unit_cell.atom_true):
                        trans_for_given_rot.append(w)
                        _same_index.append(__same_index)

            w_for_given_rot.append(trans_for_given_rot)
            same_index.append(_same_index)

        look_up_table = np.array([0, 0, 0, 0, 0, 0])  # num of following point-group operations: m, 1, 2, 3, 4, 6
        for ind_ind, _rot_ind in enumerate(rot_ind):
            if w_for_given_rot[ind_ind]:
                _W_candidate = np.identity(3)
                _W_candidate[np.ix_(self.user_arg.periodicity.nonzero()[0],
                                    self.user_arg.periodicity.nonzero()[0])] = W_candidate.copy()[_rot_ind]
                self.W_select.append(_W_candidate)
                self.w_select.append(w_for_given_rot[ind_ind])
                self.same_index_select.append(same_index[ind_ind])

                look_up = (np.trace(W_candidate[_rot_ind]), np.linalg.det(W_candidate[_rot_ind]))
                if look_up == (0.0, -1.0):
                    look_up_table[0] += 1
                elif look_up == (2.0, 1.0):
                    look_up_table[1] += 1
                elif look_up == (-2.0, 1.0):
                    look_up_table[2] += 1
                elif look_up == (-1.0, 1.0):
                    look_up_table[3] += 1
                elif look_up == (0.0, 1.0):
                    look_up_table[4] += 1
                elif look_up == (1.0, 1.0):
                    look_up_table[5] += 1
                else:
                    print('What is this operation?')
                    assert False

        if np.allclose(look_up_table, np.array([0, 1, 0, 0, 0, 0])):
            # print('Point group = 1')
            self.point_group = '1'
        elif np.allclose(look_up_table, np.array([0, 1, 1, 0, 0, 0])):
            # print('Point group = 2')
            self.point_group = '2'
        elif np.allclose(look_up_table, np.array([1, 1, 0, 0, 0, 0])):
            # print('Point group = m')
            self.point_group = 'm'
        elif np.allclose(look_up_table, np.array([2, 1, 1, 0, 0, 0])):
            # print('Point group = 2mm')
            self.point_group = '2mm'
        elif np.allclose(look_up_table, np.array([2, 2, 0, 0, 0, 0])):
            # print('Point group = m (cm)')
            self.point_group = 'm (cm)'
        elif np.allclose(look_up_table, np.array([4, 2, 2, 0, 0, 0])):
            # print('Point group = 2mm (c2mm)')
            self.point_group = '2mm (c2mm)'
        elif np.allclose(look_up_table, np.array([0, 1, 1, 0, 2, 0])):
            # print('Point group = 4')
            self.point_group = '4'
        elif np.allclose(look_up_table, np.array([4, 1, 1, 0, 2, 0])):
            # print('Point group = 4mm')
            self.point_group = '4mm'
        elif np.allclose(look_up_table, np.array([0, 1, 0, 2, 0, 0])):
            # print('Point group = 3')
            self.point_group = '3'
        elif np.allclose(look_up_table, np.array([3, 1, 0, 2, 0, 0])):
            # print('Point group = 3m')
            self.point_group = '3m'
        elif np.allclose(look_up_table, np.array([0, 1, 1, 2, 0, 2])):
            # print('Point group = 6')
            self.point_group = '6'
        elif np.allclose(look_up_table, np.array([6, 1, 1, 2, 0, 2])):
            # print('Point group = 6mm')
            self.point_group = '6mm'
        else:
            print('What is this point group?')
            assert False

        return self.W_select, self.w_select, self.same_index_select

    def search_image_atom(self):
        for _ind, _ in enumerate(self.unit_cell.atom_true):
            if self.require_atom:
                found_flag = False
                for W_ind, same in enumerate(self.same_index_select):
                    if same[0][_ind] in self.require_atom:
                        self.point_group_ind.append(W_ind)
                        self.not_require_atom.append(_ind)
                        found_flag = True
                        break
                if not found_flag:
                    self.require_atom.append(_ind)
            else:
                self.require_atom.append(_ind)

        for _W_ind, _W_select in enumerate(self.W_select):
            _same_supercell_image_index = self.search_cell_image_index(_W_select, self.super_cell)
            self.same_supercell_index_select.append(_same_supercell_image_index)

        return self.point_group_ind, self.require_atom, self.not_require_atom, self.same_supercell_index_select

    def search_cell_image_index(self, W_direct, cell):
        # conserving the image index in primitive cell
        _enlarge = int(len(cell.atom_type) / len(self.unit_cell.atom_type))

        satom_true_original = np.transpose(cell.atom_direct.copy()[cell.atom_true, :])  # [3, satom_true]
        _same_cell_index = []
        for _atom_index, _image_index in enumerate(self.same_index_select[self.find_point_group_index(W_direct)][0]):
            _satom_in_primitive = satom_true_original[:, _atom_index * _enlarge]  # [3, ]
            _satom_in_primitive_rot = satom_true_original[:, _image_index * _enlarge]  # [3, ]

            __same_cell_index = []
            for _satom_index, value in enumerate(cell.atom_true):
                _min_vector = satom_true_original[:, _satom_index] - _satom_in_primitive
                _min_vector_rot = np.dot(W_direct, _min_vector)

                for w in self.w_select[self.find_point_group_index(W_direct)]:
                    _satom_in_primitive_transform = _satom_in_primitive_rot + w.reshape([3, ])

                    same_satom_type = [ind_ for ind_, val_ in enumerate(cell.atom_true)
                                       if cell.atom_type[val_] == cell.atom_type[value]]

                    for _, same_satom_index in enumerate(same_satom_type):
                        delta_x = _satom_in_primitive_transform + _min_vector_rot - satom_true_original[:, same_satom_index]
                        delta_x_cart = np.matmul(np.transpose(cell.lattice_matrix), delta_x - np.rint(delta_x))

                        if np.allclose(delta_x_cart, np.zeros([3, ]), atol=1e-06):
                            # if same_satom_index not in __same_supercell_index:
                            __same_cell_index.append(same_satom_index)
                            #    break

                    if len(__same_cell_index) == len(cell.atom_true):
                        _same_cell_index.append(__same_cell_index)

        return _same_cell_index

    def search_self_image_atom(self):
        for _require in self.require_atom:
            _sym_point_for_require = []
            for _W_ind, _same_index_select in enumerate(self.same_index_select):
                if _same_index_select[0][_require] == _require:
                    _sym_point_for_require.append(_W_ind)
            self.point_group_for_self_require_atom.append(_sym_point_for_require)

    def search_independent_displacement(self):
        _original_basis = np.transpose(self.unit_cell.lattice_matrix.copy())
        to_cart_coord = _original_basis / np.linalg.norm(_original_basis, axis=0)
        to_direct_coord = np.linalg.inv(to_cart_coord)

        _random_direction_cart = np.array([1, 0, 1]) / np.linalg.norm(np.array([1, 0, 1]))
        for _point_group_for_self_require in self.point_group_for_self_require_atom:
            __W_displacement_cart = []
            __require_displacement_cart = []

            for __point_group_for_self_require in _point_group_for_self_require:
                tmp_W = np.identity(3)
                set_W = [tmp_W]
                tmp_W = tmp_W @ self.W_select[__point_group_for_self_require]
                while not np.allclose(tmp_W, np.identity(3)):
                    set_W.append(tmp_W)
                    tmp_W = tmp_W @ self.W_select[__point_group_for_self_require]

                set_W_in_cart = [to_cart_coord @ W_select @ to_direct_coord for W_select in set_W]

                for _W_in_cart in set_W_in_cart:
                    __W_displacement_cart.append(_W_in_cart)
                    _image_direction_cart = _W_in_cart @ _random_direction_cart.copy()
                    __require_displacement_cart.append(_image_direction_cart)

            _W_displacement_cart = [__W_displacement_cart[0]]
            _W_index = [self.find_point_group_index(to_direct_coord @ __W_displacement_cart[0] @ to_cart_coord)]
            _require_displacement_cart = [__require_displacement_cart[0]]
            for i in range(1, len(__require_displacement_cart)):
                independent = True
                for j in range(0, i):
                    # independence test by Cauchy–Schwarz inequality
                    _inner_product = np.dot(__require_displacement_cart[i], __require_displacement_cart[j])
                    _diff = np.dot(__require_displacement_cart[i], __require_displacement_cart[i]) \
                            * np.dot(__require_displacement_cart[j], __require_displacement_cart[j]) \
                            - np.dot(_inner_product, _inner_product)

                    if _diff < 1e-06:
                        independent = False
                        break

                if independent:
                    if len(_require_displacement_cart) < 3:
                        _W_displacement_cart.append(__W_displacement_cart[i])
                        _W_index.append(self.find_point_group_index(to_direct_coord @ __W_displacement_cart[i] @ to_cart_coord))
                        _require_displacement_cart.append(__require_displacement_cart[i])
                    else:
                        break

            self.independent_by_W_index.append(_W_index)
            self.independent_by_W_displacement_cart.append(_W_displacement_cart)
            self.independent_by_single_displacement_cart.append(_require_displacement_cart)

    def gen_additional_displacement(self):
        for _independent_by_single_displacement_cart in self.independent_by_single_displacement_cart:
            _independent_additional_displacement_cart = []
            additional_displacement_candidate = [np.array([0, 1, 1]) / np.linalg.norm(np.array([0, 1, 1])),
                                                 np.array([1, 1, 0]) / np.linalg.norm(np.array([1, 1, 0]))]
            if len(_independent_by_single_displacement_cart) == 3:
                pass

            elif len(_independent_by_single_displacement_cart) == 2:
                for _additional_displacement in additional_displacement_candidate:
                    independent = True
                    for __independent_by_single_displacement in _independent_by_single_displacement_cart:
                        _inner_product = np.dot(_additional_displacement, __independent_by_single_displacement)
                        _diff = np.dot(_additional_displacement, _additional_displacement) \
                                * np.dot(__independent_by_single_displacement, __independent_by_single_displacement) \
                                - np.dot(_inner_product, _inner_product)

                        if _diff < 1e-06:
                            independent = False
                            break
                    if independent:
                        _independent_additional_displacement_cart.append(_additional_displacement)
                        break

            elif len(_independent_by_single_displacement_cart) == 1:
                for _additional_displacement in additional_displacement_candidate:
                    _independent_additional_displacement_cart.append(_additional_displacement)

            self.independent_additional_displacement_cart.append(_independent_additional_displacement_cart)

    def find_point_group_index(self, W_direct):
        found_flag = False
        for _W_ind, _W_select in enumerate(self.W_select):
            if np.allclose(_W_select, W_direct, atol=1e-06):
                found_flag = True
                return _W_ind

        if not found_flag:
            assert False
