import numpy as np
from InterPhon.core import UnitCell

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

unit_cell = UnitCell()
unit_cell.read_unit_cell('./POSCAR_InAs_A')
print(unit_cell.atom_type)
print(unit_cell.atom_true)

atom_true_original = np.transpose(unit_cell.atom_direct[unit_cell.atom_true, :])

G_metric = np.dot(unit_cell.lattice_matrix[0:2, 0:2], np.transpose(unit_cell.lattice_matrix[0:2, 0:2]))  # metric tensor
# print(G_metric)

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

    w_candidate = [np.array([0.0, 0.0, 0.0])]
    # for atom_ind in range(atom_true.shape[1]):
    for index, value in enumerate(unit_cell.atom_true):
        other_atom_ind = [i for i in unit_cell.atom_true]
        other_atom_ind.remove(value)
        for other_index, other_ind in enumerate(other_atom_ind):
            if unit_cell.atom_type[other_ind] == unit_cell.atom_type[value]:
                w = np.round_(atom_true_original[:, other_index] - atom_true_rot[:, index], 6)
                if w[2] == 0.0:
                    w_, _ = np.modf(w)
                    if not np.allclose(w_, np.zeros([3, ]), atol=1e-06):
                        w_candidate.append(w_)
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
print(w_for_given_rot)
print(same_index)

look_up_table = np.array([0, 0, 0, 0, 0, 0])  # num of following point-group operations: m, 1, 2, 3, 4, 6
for ind_ind, _rot_ind in enumerate(rot_ind):
    if w_for_given_rot[ind_ind]:
        # print(W_candidate[ind])
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

not_require = []
for ind_ind, _rot_ind in enumerate(rot_ind):
    if w_for_given_rot[ind_ind]:
        for _, _same in enumerate(same_index[ind_ind]):
            for __ind, __same in enumerate(_same):
                if __ind < __same:
                    if __same not in not_require:
                        not_require.append(__same)
                        print(ind_ind)
print(not_require)
