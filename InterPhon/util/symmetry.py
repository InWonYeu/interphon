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
unit_cell.read_unit_cell('./POSCAR')
# print(unit_cell.atom_type)

atom_true = unit_cell.atom_direct[unit_cell.atom_true, :].transpose()
# print(atom_true)

G_metric = np.dot(unit_cell.lattice_matrix[0:2, 0:2], np.transpose(unit_cell.lattice_matrix[0:2, 0:2]))  # metric tensor
# print(G_metric)

rot_ind = []
for ind, rot in enumerate(W_candidate):
    G_rotate = np.dot(np.transpose(rot[0:2, 0:2]), np.dot(G_metric, rot[0:2, 0:2]))
    if np.allclose(G_metric, G_rotate, atol=1e-06):
        rot_ind.append(ind)
# print(rot_ind)

w_for_given_rot = []
for ind in rot_ind:
    # print(W_candidate[ind])
    atom_true_rot = np.dot(W_candidate[ind], atom_true)

    # print(atom_true_rot)
    # if np.allclose(atom_true_rot, atom_true):
    #     print('True')
    # else:
    #     print('False')

    w_candidate = [np.array([0.0, 0.0, 0.0])]
    for atom_ind in range(atom_true.shape[1]):
        other_atom = list(range(atom_true.shape[1]))
        other_atom.remove(atom_ind)
        for other_ind in other_atom:
            w = atom_true[:, other_ind] - atom_true_rot[:, atom_ind]
            if w[2] == 0.0:
                w_candidate.append(w)
    # print('w_candidate=', w_candidate)

    trans_for_given_rot = []
    for w in w_candidate:
        atom_transform = atom_true_rot + w.reshape([3, 1])
        delta_x = atom_transform - atom_true  # It will be replaced by atom-to-atom comparison
        delta_x_cart = np.dot(np.transpose(unit_cell.lattice_matrix), delta_x - np.rint(delta_x))

        # print(del_x.shape)
        if np.allclose(delta_x_cart, np.zeros([3, 3]), atol=1e-06):  # [3, num_atoms] It will be replaced [3, ]
            trans_for_given_rot.append(w)

    w_for_given_rot.append(trans_for_given_rot)
print(w_for_given_rot)
# print(np.array(w_for_given_rot).shape)

look_up_table = np.array([0, 0, 0, 0, 0, 0])  # num of following point-group operations: m, 1, 2, 3, 4, 6
for ind_ind, ind in enumerate(rot_ind):
    if w_for_given_rot[ind_ind]:
        print(W_candidate[ind])
        look_up = (np.trace(W_candidate[ind][0:2, 0:2]), np.linalg.det(W_candidate[ind][0:2, 0:2]))
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
    print('Point group = cm')
elif np.allclose(look_up_table, np.array([4, 2, 2, 0, 0, 0])):
    print('Point group = c2mm')
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
