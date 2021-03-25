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
print(unit_cell.atom_type)

atom_true = unit_cell.atom_direct[unit_cell.atom_true, :].transpose()
print(atom_true)

G_metric = np.dot(unit_cell.lattice_matrix[0:2, 0:2], np.transpose(unit_cell.lattice_matrix[0:2, 0:2]))  # metric tensor
print(G_metric)

rot_ind = []
for ind, rot in enumerate(W_candidate):
    G_rotate = np.dot(np.transpose(rot[0:2, 0:2]), np.dot(G_metric, rot[0:2, 0:2]))
    if np.allclose(G_metric, G_rotate, atol=1e-06):
        rot_ind.append(ind)

w_for_given_rot = []
for ind in rot_ind:
    print(W_candidate[ind])
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
        if np.allclose(delta_x_cart, np.zeros([3, 3]), atol=1e-06):  # It will be replaced by np.zeros([3, ])
            trans_for_given_rot.append(w)

    w_for_given_rot.append(trans_for_given_rot)
print(w_for_given_rot)
