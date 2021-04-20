import numpy as np
from typing import List
from InterPhon.util import KptPath
from InterPhon import error


def gamma_centered(k_file_lines: List[str], _ind_pbc) -> tuple:
    """
    Automatic generation scheme sampled by Gamma-centered grid.

    :param k_file_lines: (List[str]) List of each line of VASP KPOINTS file.
    :param _ind_pbc: (np.ndarray[int]) Indices for the direction of periodic boundary.
    :return: (tuple) A set of k_points, auto_k_points
    """
    k_points = []
    auto_k_points = [int(num_k_point) for num_k_point in k_file_lines[3].split()[0:3]]
    try:
        shift = [float(s) for s in k_file_lines[4].split()[0:3]]
    except IndexError:
        shift = [0.0, 0.0, 0.0]

    if _ind_pbc.shape[0] == 0:
        k_points.append(np.zeros((3,)))

    elif _ind_pbc.shape[0] == 1:
        __tmp_1st = np.zeros((3,))
        __tmp_1st[_ind_pbc[0]] = 1.0
        for ind_1st in range(auto_k_points[_ind_pbc[0]]):
            k_points.append(
                __tmp_1st * (ind_1st + shift[_ind_pbc[0]]) / float(auto_k_points[_ind_pbc[0]]))

            # if auto_k_points[_ind_pbc[0]] % 2 == 1:
            #     k_points.append(
            #         __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
            #         __tmp_1st * 0.5)
            # else:
            #     k_points.append(
            #         __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
            #         __tmp_1st * 0.5)

    elif _ind_pbc.shape[0] == 2:
        __tmp_1st = np.zeros((3,))
        __tmp_2nd = np.zeros((3,))
        __tmp_1st[_ind_pbc[0]] = 1.0
        __tmp_2nd[_ind_pbc[1]] = 1.0
        for ind_1st in range(auto_k_points[_ind_pbc[0]]):
            for ind_2nd in range(auto_k_points[_ind_pbc[1]]):
                k_points.append(
                    __tmp_1st * (ind_1st + shift[_ind_pbc[0]]) / float(auto_k_points[_ind_pbc[0]]) +
                    __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]]) / float(auto_k_points[_ind_pbc[1]]))

                # if auto_k_points[_ind_pbc[0]] % 2 == 1:
                #     if auto_k_points[_ind_pbc[1]] % 2 == 1:
                #         k_points.append(
                #             __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                #             __tmp_1st * 0.5 +
                #             __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                #             __tmp_2nd * 0.5)
                #     else:
                #         k_points.append(
                #             __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                #             __tmp_1st * 0.5 +
                #             __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.0) / float(auto_k_points[_ind_pbc[1]]) -
                #             __tmp_2nd * 0.5)
                # else:
                #     if auto_k_points[_ind_pbc[1]] % 2 == 1:
                #         k_points.append(
                #             __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
                #             __tmp_1st * 0.5 +
                #             __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                #             __tmp_2nd * 0.5)
                #     else:
                #         k_points.append(
                #             __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
                #             __tmp_1st * 0.5 +
                #             __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.0) / float(auto_k_points[_ind_pbc[1]]) -
                #             __tmp_2nd * 0.5)

    elif _ind_pbc.shape[0] == 3:
        __tmp_1st = np.zeros((3,))
        __tmp_2nd = np.zeros((3,))
        __tmp_3rd = np.zeros((3,))
        __tmp_1st[_ind_pbc[0]] = 1.0
        __tmp_2nd[_ind_pbc[1]] = 1.0
        __tmp_3rd[_ind_pbc[2]] = 1.0
        for ind_1st in range(auto_k_points[_ind_pbc[0]]):
            for ind_2nd in range(auto_k_points[_ind_pbc[1]]):
                for ind_3rd in range(auto_k_points[_ind_pbc[2]]):
                    k_points.append(
                        __tmp_1st * (ind_1st + shift[_ind_pbc[0]]) / float(auto_k_points[_ind_pbc[0]]) +
                        __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]]) / float(auto_k_points[_ind_pbc[1]]) +
                        __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]]) / float(auto_k_points[_ind_pbc[2]]))

                    # if auto_k_points[_ind_pbc[0]] % 2 == 1:
                    #     if auto_k_points[_ind_pbc[1]] % 2 == 1:
                    #         if auto_k_points[_ind_pbc[2]] % 2 == 1:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.5) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    #         else:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.0) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    #     else:
                    #         if auto_k_points[_ind_pbc[2]] % 2 == 1:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.0) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.5) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    #         else:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.0) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.0) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    # else:
                    #     if auto_k_points[_ind_pbc[1]] % 2 == 1:
                    #         if auto_k_points[_ind_pbc[2]] % 2 == 1:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.5) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    #         else:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.0) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    #     else:
                    #         if auto_k_points[_ind_pbc[2]] % 2 == 1:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.0) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.5) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)
                    #         else:
                    #             k_points.append(
                    #                 __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.0) / float(auto_k_points[_ind_pbc[0]]) -
                    #                 __tmp_1st * 0.5 +
                    #                 __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.0) / float(auto_k_points[_ind_pbc[1]]) -
                    #                 __tmp_2nd * 0.5 +
                    #                 __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.0) / float(auto_k_points[_ind_pbc[2]]) -
                    #                 __tmp_3rd * 0.5)

    return k_points, auto_k_points


def monkhorst_pack(k_file_lines: List[str], _ind_pbc) -> tuple:
    """
    Automatic generation scheme sampled by Monkhorst-Pack grid.

    :param k_file_lines: (List[str]) List of each line of VASP KPOINTS file.
    :param _ind_pbc: (np.ndarray[int]) Indices for the direction of periodic boundary.
    :return: (tuple) A set of k_points, auto_k_points
    """
    k_points = []
    auto_k_points = [int(num_k_point) for num_k_point in k_file_lines[3].split()[0:3]]
    try:
        shift = [float(s) for s in k_file_lines[4].split()[0:3]]
    except IndexError:
        shift = [0.0, 0.0, 0.0]

    if _ind_pbc.shape[0] == 0:
        k_points.append(np.zeros((3,)))

    elif _ind_pbc.shape[0] == 1:
        __tmp_1st = np.zeros((3,))
        __tmp_1st[_ind_pbc[0]] = 1.0
        for ind_1st in range(auto_k_points[_ind_pbc[0]]):
            k_points.append(
                __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]))

            # k_points.append(
            #     __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
            #     __tmp_1st * 0.5)

    elif _ind_pbc.shape[0] == 2:
        __tmp_1st = np.zeros((3,))
        __tmp_2nd = np.zeros((3,))
        __tmp_1st[_ind_pbc[0]] = 1.0
        __tmp_2nd[_ind_pbc[1]] = 1.0
        for ind_1st in range(auto_k_points[_ind_pbc[0]]):
            for ind_2nd in range(auto_k_points[_ind_pbc[1]]):
                k_points.append(
                    __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) +
                    __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]))

                # k_points.append(
                #     __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                #     __tmp_1st * 0.5 +
                #     __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                #     __tmp_2nd * 0.5)

    elif _ind_pbc.shape[0] == 3:
        __tmp_1st = np.zeros((3,))
        __tmp_2nd = np.zeros((3,))
        __tmp_3rd = np.zeros((3,))
        __tmp_1st[_ind_pbc[0]] = 1.0
        __tmp_2nd[_ind_pbc[1]] = 1.0
        __tmp_3rd[_ind_pbc[2]] = 1.0
        for ind_1st in range(auto_k_points[_ind_pbc[0]]):
            for ind_2nd in range(auto_k_points[_ind_pbc[1]]):
                for ind_3rd in range(auto_k_points[_ind_pbc[2]]):
                    k_points.append(
                        __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) +
                        __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) +
                        __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.5) / float(auto_k_points[_ind_pbc[2]]))

                    # k_points.append(
                    #     __tmp_1st * (ind_1st + shift[_ind_pbc[0]] + 0.5) / float(auto_k_points[_ind_pbc[0]]) -
                    #     __tmp_1st * 0.5 +
                    #     __tmp_2nd * (ind_2nd + shift[_ind_pbc[1]] + 0.5) / float(auto_k_points[_ind_pbc[1]]) -
                    #     __tmp_2nd * 0.5 +
                    #     __tmp_3rd * (ind_3rd + shift[_ind_pbc[2]] + 0.5) / float(auto_k_points[_ind_pbc[2]]) -
                    #     __tmp_3rd * 0.5)

    return k_points, auto_k_points


def line_path(k_file_lines: List[str]) -> KptPath:
    """
    Line-path generation scheme connecting high symmetry points by line-segmentation.

    :param k_file_lines: (List[str]) List of each line of VASP KPOINTS file.
    :return: (KptPath) k_points
    """
    k_points = []
    num_per_segment = int(k_file_lines[1].split()[0])
    _tmp = [k_point.split()[0:3] for k_point in k_file_lines[3:] if k_point.split()]

    if len(_tmp) % 2 != 0:
        raise error.Invalid_Line_Kpath_Error
    else:
        for ind in range(0, len(_tmp), 2):
            k_points.extend(np.linspace(np.asfarray(_tmp[ind]),
                                        np.asfarray(_tmp[ind + 1]), num_per_segment))

    return k_points


def explicit_reciprocal(k_file_lines: List[str]) -> KptPath:
    """
    Explicit list of all k-points in reciprocal coordinates.

    :param k_file_lines: (List[str]) List of each line of VASP KPOINTS file.
    :return: (KptPath) k_points
    """
    num_k_points = int(k_file_lines[1].split()[0])
    k_points = [np.asfarray(k_point.split()[0:3]) for k_point in k_file_lines[3:3+num_k_points] if k_point.split()]

    return k_points
