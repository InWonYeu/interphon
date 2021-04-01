import numpy as np
from copy import deepcopy
from InterPhon.core import UnitCell, SuperCell, PreArgument

# Searching lattice point group operations
W_candidate = [np.array([[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]]), np.array([[0, 1, 0],
                                                [-1, 0, 0],
                                                [0, 0, 1]]), np.array([[0, -1, 0],
                                                                       [1, 0, 0],
                                                                       [0, 0, 1]]), np.array([[0, -1, 0],
                                                                                              [-1, 0, 0],
                                                                                              [0, 0, 1]]),
               np.array([[0, 1, 0],
                         [1, 1, 0],
                         [0, 0, 1]]), np.array([[0, 1, 0],
                                                [1, -1, 0],
                                                [0, 0, 1]]), np.array([[0, 1, 0],
                                                                       [-1, 1, 0],
                                                                       [0, 0, 1]]), np.array([[0, -1, 0],
                                                                                              [1, 1, 0],
                                                                                              [0, 0, 1]]), np.array([[0, -1, 0],
                                                                                                                     [-1, 1, 0],
                                                                                                                     [0, 0, 1]]), np.array([[0, -1, 0],
                                                                                                                                            [1, -1, 0],
                                                                                                                                            [0, 0, 1]]), np.array([[0, 1, 0],
                                                                                                                                                                   [-1, -1, 0],
                                                                                                                                                                   [0, 0, 1]]), np.array([[0, -1, 0],
                                                                                                                                                                                          [-1, -1, 0],
                                                                                                                                                                                          [0, 0, 1]]),
               np.array([[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]]), np.array([[1, 0, 0],
                                                [0, -1, 0],
                                                [0, 0, 1]]), np.array([[-1, 0, 0],
                                                                       [0, 1, 0],
                                                                       [0, 0, 1]]), np.array([[-1, 0, 0],
                                                                                              [0, -1, 0],
                                                                                              [0, 0, 1]]),
               np.array([[1, 0, 0],
                         [1, 1, 0],
                         [0, 0, 1]]), np.array([[1, 0, 0],
                                                [1, -1, 0],
                                                [0, 0, 1]]), np.array([[1, 0, 0],
                                                                       [-1, 1, 0],
                                                                       [0, 0, 1]]), np.array([[-1, 0, 0],
                                                                                              [1, 1, 0],
                                                                                              [0, 0, 1]]), np.array([[-1, 0, 0],
                                                                                                                     [-1, 1, 0],
                                                                                                                     [0, 0, 1]]), np.array([[-1, 0, 0],
                                                                                                                                            [1, -1, 0],
                                                                                                                                            [0, 0, 1]]), np.array([[1, 0, 0],
                                                                                                                                                                   [-1, -1, 0],
                                                                                                                                                                   [0, 0, 1]]), np.array([[-1, 0, 0],
                                                                                                                                                                                          [-1, -1, 0],
                                                                                                                                                                                          [0, 0, 1]]),
               np.array([[1, 1, 0],
                         [0, 1, 0],
                         [0, 0, 1]]), np.array([[1, 1, 0],
                                                [0, -1, 0],
                                                [0, 0, 1]]), np.array([[1, -1, 0],
                                                                       [0, 1, 0],
                                                                       [0, 0, 1]]), np.array([[-1, 1, 0],
                                                                                              [0, 1, 0],
                                                                                              [0, 0, 1]]), np.array([[-1, -1, 0],
                                                                                                                     [0, 1, 0],
                                                                                                                     [0, 0, 1]]), np.array([[-1, 1, 0],
                                                                                                                                            [0, -1, 0],
                                                                                                                                            [0, 0, 1]]), np.array([[1, -1, 0],
                                                                                                                                                                   [0, -1, 0],
                                                                                                                                                                   [0, 0, 1]]), np.array([[-1, -1, 0],
                                                                                                                                                                                          [0, -1, 0],
                                                                                                                                                                                          [0, 0, 1]]),
               np.array([[1, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]]), np.array([[1, 1, 0],
                                                [-1, 0, 0],
                                                [0, 0, 1]]), np.array([[1, -1, 0],
                                                                       [1, 0, 0],
                                                                       [0, 0, 1]]), np.array([[-1, 1, 0],
                                                                                              [1, 0, 0],
                                                                                              [0, 0, 1]]), np.array([[-1, -1, 0],
                                                                                                                     [1, 0, 0],
                                                                                                                     [0, 0, 1]]), np.array([[-1, 1, 0],
                                                                                                                                            [-1, 0, 0],
                                                                                                                                            [0, 0, 1]]), np.array([[1, -1, 0],
                                                                                                                                                                   [-1, 0, 0],
                                                                                                                                                                   [0, 0, 1]]), np.array([[-1, -1, 0],
                                                                                                                                                                                          [-1, 0, 0],
                                                                                                                                                                                          [0, 0, 1]]),
               ]


def symmetry_2d(unit_cell=UnitCell(), super_cell=SuperCell(), user_arg=PreArgument()):
    _enlarge = 1
    for ind, value in enumerate(user_arg.periodicity):
        if value:
            _enlarge = _enlarge * user_arg.enlargement[ind]

    _super_direct = np.empty([_enlarge, 3])
    k = 0
    for x in range(user_arg.enlargement[0]):
        for y in range(user_arg.enlargement[1]):
            for z in range(user_arg.enlargement[2]):
                _super_direct[k, 0:3] = np.array([float(x) / user_arg.enlargement[0],
                                                  float(y) / user_arg.enlargement[1],
                                                  float(z) / user_arg.enlargement[2]])
                k = k + 1
    __super_direct = _super_direct[0, 0:3]

    # metric tensor
    atom_true_original = np.transpose(unit_cell.atom_direct[unit_cell.atom_true, :])
    G_metric = np.dot(unit_cell.lattice_matrix[0:2, 0:2], np.transpose(unit_cell.lattice_matrix[0:2, 0:2]))

    rot_ind = []
    for ind, rot in enumerate(W_candidate):
        G_rotate = np.dot(np.transpose(rot[0:2, 0:2]), np.dot(G_metric, rot[0:2, 0:2]))
        if np.allclose(G_metric, G_rotate, atol=1e-06):
            rot_ind.append(ind)
    # print(rot_ind)

    w_for_given_rot = []
    same_index = []
    for ind in rot_ind:
        atom_true_rot = np.dot(W_candidate[ind], atom_true_original)

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
            for index, value in enumerate(unit_cell.atom_true):
                same_atom_type = [ind_ for ind_, val_ in enumerate(unit_cell.atom_true)
                                  if unit_cell.atom_type[val_] == unit_cell.atom_type[value]]

                for _, same_atom_index in enumerate(same_atom_type):
                    delta_x = atom_transform[:, index] - atom_true_original[:, same_atom_index]  # atom-to-atom comparison
                    delta_x_cart = np.matmul(np.transpose(unit_cell.lattice_matrix), delta_x - np.rint(delta_x))

                    if np.allclose(delta_x_cart, np.zeros([3, ]), atol=1e-06):
                        if same_atom_index not in __same_index:
                            __same_index.append(same_atom_index)
                            break

                if len(__same_index) == len(unit_cell.atom_true):
                    trans_for_given_rot.append(w)
                    _same_index.append(__same_index)

        w_for_given_rot.append(trans_for_given_rot)
        same_index.append(_same_index)

    W_select = []
    w_select = []
    same_index_select = []
    look_up_table = np.array([0, 0, 0, 0, 0, 0])  # num of following point-group operations: m, 1, 2, 3, 4, 6
    for ind_ind, _rot_ind in enumerate(rot_ind):
        if w_for_given_rot[ind_ind]:
            W_select.append(W_candidate[_rot_ind])
            w_select.append(w_for_given_rot[ind_ind])
            same_index_select.append(same_index[ind_ind])

            look_up = (np.trace(W_candidate[_rot_ind][0:2, 0:2]), np.linalg.det(W_candidate[_rot_ind][0:2, 0:2]))
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
        print('Point group = 1')
    elif np.allclose(look_up_table, np.array([0, 1, 1, 0, 0, 0])):
        print('Point group = 2')
    elif np.allclose(look_up_table, np.array([1, 1, 0, 0, 0, 0])):
        print('Point group = m')
    elif np.allclose(look_up_table, np.array([2, 1, 1, 0, 0, 0])):
        print('Point group = 2mm')
    elif np.allclose(look_up_table, np.array([2, 2, 0, 0, 0, 0])):
        print('Point group = m (cm)')
    elif np.allclose(look_up_table, np.array([4, 2, 2, 0, 0, 0])):
        print('Point group = 2mm (c2mm)')
    elif np.allclose(look_up_table, np.array([0, 1, 1, 0, 2, 0])):
        print('Point group = 4')
    elif np.allclose(look_up_table, np.array([4, 1, 1, 0, 2, 0])):
        print('Point group = 4mm')
    elif np.allclose(look_up_table, np.array([0, 1, 0, 2, 0, 0])):
        print('Point group = 3')
    elif np.allclose(look_up_table, np.array([3, 1, 0, 2, 0, 0])):
        print('Point group = 3m')
    elif np.allclose(look_up_table, np.array([0, 1, 1, 2, 0, 2])):
        print('Point group = 6')
    elif np.allclose(look_up_table, np.array([6, 1, 1, 2, 0, 2])):
        print('Point group = 6mm')
    else:
        print('What is this point group?')
        assert False

    require = []
    not_require = []
    point_group_ind = []

    for _ind, _ in enumerate(unit_cell.atom_true):
        if require:
            found_flag = False
            for W_ind, same in enumerate(same_index_select):  # 여기서 same_supercell_index_select 도 같이 shuffle
                if same[0][_ind] in require:
                    point_group_ind.append(W_ind)
                    not_require.append(_ind)
                    found_flag = True
                    break
            if not found_flag:
                require.append(_ind)
        else:
            require.append(_ind)

    same_supercell_index_select = []
    for _W_ind, _W_select in enumerate(W_select):
        _same_supercell_index = []

        # find nearest-neighbor for each atom in the primitive unit cell
        for _atom_index, _atom_true in enumerate(unit_cell.atom_true):
            satom_true_original = np.transpose(super_cell.atom_direct[super_cell.atom_true, :]).copy()  # [3, satom_true]

            _atom_in_primitive = unit_cell.atom_direct[_atom_true] / user_arg.enlargement  # [3, ] matrix
            __same_supercell_index = []
            for _satom_index, value in enumerate(super_cell.atom_true):
                _distance = satom_true_original[:, _satom_index] - _atom_in_primitive
                for _first in (-1.0, 0.0, 1.0):
                    for _second in (-1.0, 0.0, 1.0):
                        _tmp_distance = _distance + _first * np.array([1, 0, 0]) + _second * np.array([0, 1, 0])
                        if np.dot(_distance, _distance) > np.dot(_tmp_distance, _tmp_distance):
                            satom_true_original[:, _satom_index] = satom_true_original[:, _satom_index].copy() \
                                                                   + _first * np.array([1, 0, 0]) \
                                                                   + _second * np.array([0, 1, 0])
                _min_vector = satom_true_original[:, _satom_index] - _atom_in_primitive

                _atom_in_primitive_rot = np.dot(_W_select, _atom_in_primitive)
                _min_vector_rot = np.dot(_W_select, _min_vector)
                for w in w_select[_W_ind]:
                    _atom_in_primitive_transform = _atom_in_primitive_rot + w.reshape([3, 1])

                    same_satom_type = [ind_ for ind_, val_ in enumerate(super_cell.atom_true)
                                       if super_cell.atom_type[val_] == super_cell.atom_type[value]]

                    _num_found = 0
                    for _, same_satom_index in enumerate(same_satom_type):
                        delta_x = _atom_in_primitive_rot + _min_vector_rot - satom_true_original[:, same_satom_index]
                        delta_x_cart = np.matmul(np.transpose(super_cell.lattice_matrix), delta_x - np.rint(delta_x))

                        if np.allclose(delta_x_cart, np.zeros([3, ]), atol=1e-06):
                            if same_satom_index not in __same_supercell_index:
                                __same_supercell_index.append(same_satom_index)
                                break

                    if len(__same_supercell_index) == len(super_cell.atom_true):
                        _same_supercell_index.append(__same_supercell_index)

        previous = _same_supercell_index[0]
        for _test in _same_supercell_index:
            if not np.allclose(np.array(previous), np.array(_test)):
                assert False
            previous = _test
        print("Done without error")
        same_supercell_index_select.append(_same_supercell_index)

    return W_select, w_select, same_index_select, point_group_ind, require, not_require, same_supercell_index_select
