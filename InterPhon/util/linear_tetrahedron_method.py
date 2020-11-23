"""
The demonstrated formula below is described in Supplementary Information.
For more details, see the following references:

Ref 1) Extensions of the tetrahedron method for evaluating spectral properties of solids, 
Journal of Physics C 12, 2991 (1979).
Ref 2) FermiSurfer: Fermi-surface viewer providing multiple representation schemes,
Computer Physics Communications 239, 197 (2019)
"""
import numpy as np


def tetrahedron_1d(_freq, _pdos, k_points, w_q, v_q, num_k_points, _ind_pbc):
    """
    1D version of linear tetrahedron method: the linear line method.
    The detailed demonstration of the formula is written in Supplementary Information.

    :param _freq: (np.ndarray[float]) '(num_dos,) size' frequencies where DOS will be evaluated.
    :param _pdos: (np.ndarray[float]) '(num_mode, num_dos) size' partial DOS for each phonon mode.
    :param k_points: (List[np.ndarray]) List of '(3,) size' k-point.
    :param w_q: (np.ndarray[float]) '(num_k_points, num_mode) size' eigenfrequency.
    :param v_q: (np.ndarray[float]) '(num_k_points, num_mode, num_mode) size' eigenmode.
    :param num_k_points: (List[int]) '(3,) size' list for the number of automatic grid.
    :param _ind_pbc: (np.ndarray[int]) Indices for the direction of periodic boundary.
    :return: _pdos
    """
    #  1st Brillouin Zone is segmented by line components connecting the sampled k-points grid.
    #  Integration is performed for each segmentation by linear interpolation.
    #
    #  line components:
    #  0-------–1
    #
    #  line 1: 0-1
    num_1st_kpt = num_k_points[_ind_pbc[0]]
    length_bz = 1.0 ** 1
    length_ith = np.sqrt(np.dot(k_points[1 % num_1st_kpt] - k_points[0], k_points[1 % num_1st_kpt] - k_points[0]))

    for ind_x, __freq in enumerate(_freq):
        for ind_freq in range(v_q.shape[1]):
            _w_q = w_q[:, ind_freq]
            _v_q = v_q[:, ind_freq, :]
            for ind_1st_kpt in range(num_1st_kpt):
                # line 1: 0-1
                # for ind_mode in range(v_q.shape[2]):
                __w0 = _w_q[ind_1st_kpt]
                __F0 = abs(_v_q[ind_1st_kpt, :]) ** 2
                __w1 = _w_q[(ind_1st_kpt + 1) % num_1st_kpt]
                __F1 = abs(_v_q[(ind_1st_kpt + 1) % num_1st_kpt, :]) ** 2

                __w0, __w1 = np.asfarray([__w0, __w1])[np.argsort(np.asfarray([__w0, __w1]))]
                __F0, __F1 = np.asfarray([__F0, __F1])[np.argsort(np.asfarray([__w0, __w1])), :]

                if __w0 <= __freq < __w1:
                    __g = 1.0 / (__w1 - __w0)

                    __I0 = (__freq - __w1) / (__w0 - __w1)
                    __I1 = (__freq - __w0) / (__w1 - __w0)

                    _pdos[:, ind_x] += length_ith / length_bz * __g \
                                                 * (__I0 * __F0 + __I1 * __F1)

                else:
                    pass

    return _pdos


def tetrahedron_2d(_freq, _pdos, k_points, w_q, v_q, num_k_points, _ind_pbc):
    """
    2D version of linear tetrahedron method: the linear triangle method.
    The detailed demonstration of the formula is written in Supplementary Information.

    :param _freq: (np.ndarray[float]) '(num_dos,) size' frequencies where DOS will be evaluated.
    :param _pdos: (np.ndarray[float]) '(num_mode, num_dos) size' partial DOS for each phonon mode.
    :param k_points: (List[np.ndarray]) List of '(3,) size' k-point.
    :param w_q: (np.ndarray[float]) '(num_k_points, num_mode) size' eigenfrequency.
    :param v_q: (np.ndarray[float]) '(num_k_points, num_mode, num_mode) size' eigenmode.
    :param num_k_points: (List[int]) '(3,) size' list for the number of automatic grid.
    :param _ind_pbc: (np.ndarray[int]) Indices for the direction of periodic boundary.
    :return: _pdos
    """
    #  1st Brillouin Zone is segmented by triangle components connecting the sampled k-points grid.
    #  Integration is performed for each segmentation by linear interpolation.
    #
    #  triangle components:
    #  0--------1
    #  |        |
    #  |        |
    #  |        |
    #  2-------–3
    #
    # triangle 1: 0-1-3
    # triangle 2: 0-2-3
    num_1st_kpt = num_k_points[_ind_pbc[0]]
    num_2nd_kpt = num_k_points[_ind_pbc[1]]

    __k01 = k_points[1 % num_2nd_kpt] - k_points[0]
    __k02 = k_points[num_2nd_kpt] - k_points[0]

    __tmp = np.ones((3,))
    __tmp[_ind_pbc[0]] = 0.0
    __tmp[_ind_pbc[1]] = 0.0

    area_bz = 1.0 ** 2
    area_ith = 1.0 / 2.0 * abs(np.dot(__tmp, np.cross(__k01, __k02)))

    for ind_x, __freq in enumerate(_freq):
        for ind_freq in range(v_q.shape[1]):
            _w_q = w_q[:, ind_freq]
            _v_q = v_q[:, ind_freq, :]
            for ind_1st_kpt in range(num_1st_kpt):
                for ind_2nd_kpt in range(num_2nd_kpt):
                    # triangle 1: 0-1-3
                    # for ind_mode in range(v_q.shape[2]):
                    __w0 = _w_q[ind_1st_kpt * num_2nd_kpt +
                                ind_2nd_kpt]
                    __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt +
                                    ind_2nd_kpt, :]) ** 2
                    __w1 = _w_q[ind_1st_kpt * num_2nd_kpt +
                                (ind_2nd_kpt + 1) % num_2nd_kpt]
                    __F1 = abs(_v_q[ind_1st_kpt * num_2nd_kpt +
                                    (ind_2nd_kpt + 1) % num_2nd_kpt, :]) ** 2
                    __w2 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt +
                                (ind_2nd_kpt + 1) % num_2nd_kpt]
                    __F2 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt +
                                    (ind_2nd_kpt + 1) % num_2nd_kpt, :]) ** 2

                    __w0, __w1, __w2 = np.asfarray([__w0, __w1, __w2])[np.argsort(np.asfarray([__w0, __w1, __w2]))]
                    __F0, __F1, __F2 = np.asfarray([__F0, __F1, __F2])[np.argsort(np.asfarray([__w0, __w1, __w2])), :]

                    if __w0 <= __freq < __w1:
                        __g = 2 * (__freq - __w0) / (__w1 - __w0) / (__w2 - __w0)

                        __I0 = 1.0 / 2.0 * ((__freq - __w1) / (__w0 - __w1) + (__freq - __w2) / (__w0 - __w2))
                        __I1 = 1.0 / 2.0 * (__freq - __w0) / (__w1 - __w0)
                        __I2 = 1.0 / 2.0 * (__freq - __w0) / (__w2 - __w0)

                        _pdos[:, ind_x] += area_ith / area_bz * __g \
                                                     * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2)

                    elif __w1 <= __freq < __w2:
                        __g = 2 * (__w2 - __freq) / (__w2 - __w1) / (__w2 - __w0)

                        __I0 = 1.0 / 2.0 * (__freq - __w2) / (__w0 - __w2)
                        __I1 = 1.0 / 2.0 * (__freq - __w2) / (__w1 - __w2)
                        __I2 = 1.0 / 2.0 * ((__freq - __w0) / (__w2 - __w0) + (__freq - __w1) / (__w2 - __w1))

                        _pdos[:, ind_x] += area_ith / area_bz * __g \
                                                     * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2)

                    else:
                        pass

                    # triangle 2: 0-2-3
                    # for ind_mode in range(v_q.shape[2]):
                    __w0 = _w_q[ind_1st_kpt * num_2nd_kpt +
                                ind_2nd_kpt]
                    __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt +
                                    ind_2nd_kpt, :]) ** 2
                    __w1 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt +
                                ind_2nd_kpt]
                    __F1 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt +
                                    ind_2nd_kpt, :]) ** 2
                    __w2 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt +
                                (ind_2nd_kpt + 1) % num_2nd_kpt]
                    __F2 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt +
                                    (ind_2nd_kpt + 1) % num_2nd_kpt, :]) ** 2

                    __w0, __w1, __w2 = np.asfarray([__w0, __w1, __w2])[np.argsort(np.asfarray([__w0, __w1, __w2]))]
                    __F0, __F1, __F2 = np.asfarray([__F0, __F1, __F2])[np.argsort(np.asfarray([__w0, __w1, __w2])), :]

                    if __w0 <= __freq < __w1:
                        __g = 2 * (__freq - __w0) / (__w1 - __w0) / (__w2 - __w0)

                        __I0 = 1.0 / 2.0 * ((__freq - __w1) / (__w0 - __w1) + (__freq - __w2) / (__w0 - __w2))
                        __I1 = 1.0 / 2.0 * (__freq - __w0) / (__w1 - __w0)
                        __I2 = 1.0 / 2.0 * (__freq - __w0) / (__w2 - __w0)

                        _pdos[:, ind_x] += area_ith / area_bz * __g \
                                                     * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2)

                    elif __w1 <= __freq < __w2:
                        __g = 2 * (__w2 - __freq) / (__w2 - __w1) / (__w2 - __w0)

                        __I0 = 1.0 / 2.0 * (__freq - __w2) / (__w0 - __w2)
                        __I1 = 1.0 / 2.0 * (__freq - __w2) / (__w1 - __w2)
                        __I2 = 1.0 / 2.0 * ((__freq - __w0) / (__w2 - __w0) + (__freq - __w1) / (__w2 - __w1))

                        _pdos[:, ind_x] += area_ith / area_bz * __g \
                                                     * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2)

                    else:
                        pass

    return _pdos


def tetrahedron_3d(_freq, _pdos, k_points, w_q, v_q, num_k_points, _ind_pbc):
    """
    3D version of linear tetrahedron method: the linear tetrahedron method.
    The detailed demonstration of the formula is written in Supplementary Information.

    :param _freq: (np.ndarray[float]) '(num_dos,) size' frequencies where DOS will be evaluated.
    :param _pdos: (np.ndarray[float]) '(num_mode, num_dos) size' partial DOS for each phonon mode.
    :param k_points: (List[np.ndarray]) List of '(3,) size' k-point.
    :param w_q: (np.ndarray[float]) '(num_k_points, num_mode) size' eigenfrequency.
    :param v_q: (np.ndarray[float]) '(num_k_points, num_mode, num_mode) size' eigenmode.
    :param num_k_points: (List[int]) '(3,) size' list for the number of automatic grid.
    :param _ind_pbc: (np.ndarray[int]) Indices for the direction of periodic boundary.
    :return: _pdos
    """
    #  1st Brillouin Zone is segmented by tetrahedron components connecting the sampled k-points grid.
    #  Integration is performed for each segmentation by linear interpolation.
    #
    #  tetrahedron components:
    #     4--------5
    #    /|       /|
    #   / |      / |
    #  6--------7  |
    #  |  0-----|--1
    #  | /      | /
    #  |/       |/
    #  2--------3
    #
    # tetrahedron 1: 0-1-5-7
    # tetrahedron 2: 0-1-3-7
    # tetrahedron 3: 0-2-3-7
    # tetrahedron 4: 0-4-5-7
    # tetrahedron 5: 0-4-6-7
    # tetrahedron 6: 0-2-6-7
    num_1st_kpt = num_k_points[_ind_pbc[0]]
    num_2nd_kpt = num_k_points[_ind_pbc[1]]
    num_3rd_kpt = num_k_points[_ind_pbc[2]]

    __k01 = k_points[1 % num_2nd_kpt] - k_points[0]
    __k02 = k_points[num_2nd_kpt] - k_points[0]
    __k03 = k_points[num_2nd_kpt * num_3rd_kpt] - k_points[0]

    volume_bz = 1.0 ** 3
    volume_ith = 1.0 / 6.0 * abs(np.dot(__k01, np.cross(__k02, __k03)))

    for ind_x, __freq in enumerate(_freq):
        for ind_freq in range(v_q.shape[1]):
            _w_q = w_q[:, ind_freq]
            _v_q = v_q[:, ind_freq, :]
            for ind_1st_kpt in range(num_1st_kpt):
                for ind_2nd_kpt in range(num_2nd_kpt):
                    for ind_3rd_kpt in range(num_3rd_kpt):
                        # tetrahedron 1: 0-1-5-7
                        # for ind_mode in range(v_q.shape[2]):
                        __w0 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w1 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F1 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2
                        __w2 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F2 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2
                        __w3 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F3 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2

                        __w0, __w1, __w2, __w3 = np.asfarray([__w0,
                                                              __w1,
                                                              __w2,
                                                              __w3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3]))]
                        __F0, __F1, __F2, __F3 = np.asfarray([__F0,
                                                              __F1,
                                                              __F2,
                                                              __F3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3])), :]

                        if __w0 <= __freq < __w1:
                            __g = 3 * (__freq - __w0) ** 2 / (__w1 - __w0) / (__w2 - __w0) / (__w3 - __w0)

                            __I0 = 1.0 / 3.0 * ((__freq - __w1) / (__w0 - __w1)
                                                + (__freq - __w2) / (__w0 - __w2)
                                                + (__freq - __w3) / (__w0 - __w3))
                            __I1 = 1.0 / 3.0 * (__freq - __w0) / (__w1 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w0) / (__w2 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                                         * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w1 <= __freq < __w2:
                            __g = 3 / (__w3 - __w0) \
                                  * ((__freq - __w1) * (__freq - __w3) / ((__w2 - __w1) * (__w1 - __w3))
                                     + (__freq - __w0) * (__freq - __w2) / ((__w2 - __w0) * (__w1 - __w2)))

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3) \
                                   + (__freq - __w2) / (__w0 - __w2) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I1 = 1.0 / 3.0 * (__freq - __w2) / (__w1 - __w2) \
                                   + (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w1) / (__w2 - __w1) \
                                   + (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0) \
                                   + (__freq - __w1) / (__w3 - __w1) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                                         * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w2 <= __freq < __w3:
                            __g = 3 * (__w3 - __freq) ** 2 / (__w3 - __w0) / (__w3 - __w1) / (__w3 - __w2)

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3)
                            __I1 = 1.0 / 3.0 * (__freq - __w3) / (__w1 - __w3)
                            __I2 = 1.0 / 3.0 * (__freq - __w3) / (__w2 - __w3)
                            __I3 = 1.0 / 3.0 * ((__freq - __w0) / (__w3 - __w0)
                                                + (__freq - __w1) / (__w3 - __w1)
                                                + (__freq - __w2) / (__w3 - __w2))

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                                         * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        else:
                            pass

                        # tetrahedron 2: 0-1-3-7
                        # for ind_mode in range(v_q.shape[2]):
                        __w0 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w1 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F1 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2
                        __w2 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F2 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2
                        __w3 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F3 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2

                        __w0, __w1, __w2, __w3 = np.asfarray([__w0,
                                                              __w1,
                                                              __w2,
                                                              __w3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3]))]
                        __F0, __F1, __F2, __F3 = np.asfarray([__F0,
                                                              __F1,
                                                              __F2,
                                                              __F3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3])), :]

                        if __w0 <= __freq < __w1:
                            __g = 3 * (__freq - __w0) ** 2 / (__w1 - __w0) / (__w2 - __w0) / (__w3 - __w0)

                            __I0 = 1.0 / 3.0 * ((__freq - __w1) / (__w0 - __w1)
                                                + (__freq - __w2) / (__w0 - __w2)
                                                + (__freq - __w3) / (__w0 - __w3))
                            __I1 = 1.0 / 3.0 * (__freq - __w0) / (__w1 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w0) / (__w2 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w1 <= __freq < __w2:
                            __g = 3 / (__w3 - __w0) \
                                  * ((__freq - __w1) * (__freq - __w3) / ((__w2 - __w1) * (__w1 - __w3))
                                     + (__freq - __w0) * (__freq - __w2) / ((__w2 - __w0) * (__w1 - __w2)))

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3) \
                                   + (__freq - __w2) / (__w0 - __w2) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I1 = 1.0 / 3.0 * (__freq - __w2) / (__w1 - __w2) \
                                   + (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w1) / (__w2 - __w1) \
                                   + (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0) \
                                   + (__freq - __w1) / (__w3 - __w1) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w2 <= __freq < __w3:
                            __g = 3 * (__w3 - __freq) ** 2 / (__w3 - __w0) / (__w3 - __w1) / (__w3 - __w2)

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3)
                            __I1 = 1.0 / 3.0 * (__freq - __w3) / (__w1 - __w3)
                            __I2 = 1.0 / 3.0 * (__freq - __w3) / (__w2 - __w3)
                            __I3 = 1.0 / 3.0 * ((__freq - __w0) / (__w3 - __w0)
                                                + (__freq - __w1) / (__w3 - __w1)
                                                + (__freq - __w2) / (__w3 - __w2))

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        else:
                            pass

                        # tetrahedron 3: 0-2-3-7
                        # for ind_mode in range(v_q.shape[2]):
                        __w0 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w1 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F1 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w2 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F2 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2
                        __w3 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F3 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2

                        __w0, __w1, __w2, __w3 = np.asfarray([__w0,
                                                              __w1,
                                                              __w2,
                                                              __w3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3]))]
                        __F0, __F1, __F2, __F3 = np.asfarray([__F0,
                                                              __F1,
                                                              __F2,
                                                              __F3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3])), :]

                        if __w0 <= __freq < __w1:
                            __g = 3 * (__freq - __w0) ** 2 / (__w1 - __w0) / (__w2 - __w0) / (__w3 - __w0)

                            __I0 = 1.0 / 3.0 * ((__freq - __w1) / (__w0 - __w1)
                                                + (__freq - __w2) / (__w0 - __w2)
                                                + (__freq - __w3) / (__w0 - __w3))
                            __I1 = 1.0 / 3.0 * (__freq - __w0) / (__w1 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w0) / (__w2 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w1 <= __freq < __w2:
                            __g = 3 / (__w3 - __w0) \
                                  * ((__freq - __w1) * (__freq - __w3) / ((__w2 - __w1) * (__w1 - __w3))
                                     + (__freq - __w0) * (__freq - __w2) / ((__w2 - __w0) * (__w1 - __w2)))

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3) \
                                   + (__freq - __w2) / (__w0 - __w2) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I1 = 1.0 / 3.0 * (__freq - __w2) / (__w1 - __w2) \
                                   + (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w1) / (__w2 - __w1) \
                                   + (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0) \
                                   + (__freq - __w1) / (__w3 - __w1) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w2 <= __freq < __w3:
                            __g = 3 * (__w3 - __freq) ** 2 / (__w3 - __w0) / (__w3 - __w1) / (__w3 - __w2)

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3)
                            __I1 = 1.0 / 3.0 * (__freq - __w3) / (__w1 - __w3)
                            __I2 = 1.0 / 3.0 * (__freq - __w3) / (__w2 - __w3)
                            __I3 = 1.0 / 3.0 * ((__freq - __w0) / (__w3 - __w0)
                                                + (__freq - __w1) / (__w3 - __w1)
                                                + (__freq - __w2) / (__w3 - __w2))

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        else:
                            pass

                        # tetrahedron 4: 0-4-5-7
                        # for ind_mode in range(v_q.shape[2]):
                        __w0 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w1 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F1 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w2 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F2 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2
                        __w3 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F3 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2

                        __w0, __w1, __w2, __w3 = np.asfarray([__w0,
                                                              __w1,
                                                              __w2,
                                                              __w3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3]))]
                        __F0, __F1, __F2, __F3 = np.asfarray([__F0,
                                                              __F1,
                                                              __F2,
                                                              __F3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3])), :]

                        if __w0 <= __freq < __w1:
                            __g = 3 * (__freq - __w0) ** 2 / (__w1 - __w0) / (__w2 - __w0) / (__w3 - __w0)

                            __I0 = 1.0 / 3.0 * ((__freq - __w1) / (__w0 - __w1)
                                                + (__freq - __w2) / (__w0 - __w2)
                                                + (__freq - __w3) / (__w0 - __w3))
                            __I1 = 1.0 / 3.0 * (__freq - __w0) / (__w1 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w0) / (__w2 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w1 <= __freq < __w2:
                            __g = 3 / (__w3 - __w0) \
                                  * ((__freq - __w1) * (__freq - __w3) / ((__w2 - __w1) * (__w1 - __w3))
                                     + (__freq - __w0) * (__freq - __w2) / ((__w2 - __w0) * (__w1 - __w2)))

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3) \
                                   + (__freq - __w2) / (__w0 - __w2) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I1 = 1.0 / 3.0 * (__freq - __w2) / (__w1 - __w2) \
                                   + (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w1) / (__w2 - __w1) \
                                   + (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0) \
                                   + (__freq - __w1) / (__w3 - __w1) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w2 <= __freq < __w3:
                            __g = 3 * (__w3 - __freq) ** 2 / (__w3 - __w0) / (__w3 - __w1) / (__w3 - __w2)

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3)
                            __I1 = 1.0 / 3.0 * (__freq - __w3) / (__w1 - __w3)
                            __I2 = 1.0 / 3.0 * (__freq - __w3) / (__w2 - __w3)
                            __I3 = 1.0 / 3.0 * ((__freq - __w0) / (__w3 - __w0)
                                                + (__freq - __w1) / (__w3 - __w1)
                                                + (__freq - __w2) / (__w3 - __w2))

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        else:
                            pass

                        # tetrahedron 5: 0-4-6-7
                        # for ind_mode in range(v_q.shape[2]):
                        __w0 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w1 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F1 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w2 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F2 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w3 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F3 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2

                        __w0, __w1, __w2, __w3 = np.asfarray([__w0,
                                                              __w1,
                                                              __w2,
                                                              __w3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3]))]
                        __F0, __F1, __F2, __F3 = np.asfarray([__F0,
                                                              __F1,
                                                              __F2,
                                                              __F3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3])), :]

                        if __w0 <= __freq < __w1:
                            __g = 3 * (__freq - __w0) ** 2 / (__w1 - __w0) / (__w2 - __w0) / (__w3 - __w0)

                            __I0 = 1.0 / 3.0 * ((__freq - __w1) / (__w0 - __w1)
                                                + (__freq - __w2) / (__w0 - __w2)
                                                + (__freq - __w3) / (__w0 - __w3))
                            __I1 = 1.0 / 3.0 * (__freq - __w0) / (__w1 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w0) / (__w2 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w1 <= __freq < __w2:
                            __g = 3 / (__w3 - __w0) \
                                  * ((__freq - __w1) * (__freq - __w3) / ((__w2 - __w1) * (__w1 - __w3))
                                     + (__freq - __w0) * (__freq - __w2) / ((__w2 - __w0) * (__w1 - __w2)))

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3) \
                                   + (__freq - __w2) / (__w0 - __w2) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I1 = 1.0 / 3.0 * (__freq - __w2) / (__w1 - __w2) \
                                   + (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w1) / (__w2 - __w1) \
                                   + (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0) \
                                   + (__freq - __w1) / (__w3 - __w1) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w2 <= __freq < __w3:
                            __g = 3 * (__w3 - __freq) ** 2 / (__w3 - __w0) / (__w3 - __w1) / (__w3 - __w2)

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3)
                            __I1 = 1.0 / 3.0 * (__freq - __w3) / (__w1 - __w3)
                            __I2 = 1.0 / 3.0 * (__freq - __w3) / (__w2 - __w3)
                            __I3 = 1.0 / 3.0 * ((__freq - __w0) / (__w3 - __w0)
                                                + (__freq - __w1) / (__w3 - __w1)
                                                + (__freq - __w2) / (__w3 - __w2))

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        else:
                            pass

                        # tetrahedron 6: 0-2-6-7
                        # for ind_mode in range(v_q.shape[2]):
                        __w0 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ind_2nd_kpt * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F0 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ind_2nd_kpt * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w1 = _w_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F1 = abs(_v_q[ind_1st_kpt * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w2 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    ind_3rd_kpt]
                        __F2 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        ind_3rd_kpt, :]) ** 2
                        __w3 = _w_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                    ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                    (ind_3rd_kpt + 1) % num_3rd_kpt]
                        __F3 = abs(_v_q[((ind_1st_kpt + 1) % num_1st_kpt) * num_2nd_kpt * num_3rd_kpt +
                                        ((ind_2nd_kpt + 1) % num_2nd_kpt) * num_3rd_kpt +
                                        (ind_3rd_kpt + 1) % num_3rd_kpt, :]) ** 2

                        __w0, __w1, __w2, __w3 = np.asfarray([__w0,
                                                              __w1,
                                                              __w2,
                                                              __w3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3]))]
                        __F0, __F1, __F2, __F3 = np.asfarray([__F0,
                                                              __F1,
                                                              __F2,
                                                              __F3])[np.argsort(np.asfarray([__w0,
                                                                                             __w1,
                                                                                             __w2,
                                                                                             __w3])), :]

                        if __w0 <= __freq < __w1:
                            __g = 3 * (__freq - __w0) ** 2 / (__w1 - __w0) / (__w2 - __w0) / (__w3 - __w0)

                            __I0 = 1.0 / 3.0 * ((__freq - __w1) / (__w0 - __w1)
                                                + (__freq - __w2) / (__w0 - __w2)
                                                + (__freq - __w3) / (__w0 - __w3))
                            __I1 = 1.0 / 3.0 * (__freq - __w0) / (__w1 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w0) / (__w2 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w1 <= __freq < __w2:
                            __g = 3 / (__w3 - __w0) \
                                  * ((__freq - __w1) * (__freq - __w3) / ((__w2 - __w1) * (__w1 - __w3))
                                     + (__freq - __w0) * (__freq - __w2) / ((__w2 - __w0) * (__w1 - __w2)))

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3) \
                                   + (__freq - __w2) / (__w0 - __w2) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I1 = 1.0 / 3.0 * (__freq - __w2) / (__w1 - __w2) \
                                   + (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)
                            __I2 = 1.0 / 3.0 * (__freq - __w1) / (__w2 - __w1) \
                                   + (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w0) / (__w2 - __w0) \
                                   * (__freq - __w2) / (__w1 - __w2) / __g / (__w3 - __w0)
                            __I3 = 1.0 / 3.0 * (__freq - __w0) / (__w3 - __w0) \
                                   + (__freq - __w1) / (__w3 - __w1) \
                                   * (__freq - __w3) / (__w1 - __w3) \
                                   * (__freq - __w1) / (__w2 - __w1) / __g / (__w3 - __w0)

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        elif __w2 <= __freq < __w3:
                            __g = 3 * (__w3 - __freq) ** 2 / (__w3 - __w0) / (__w3 - __w1) / (__w3 - __w2)

                            __I0 = 1.0 / 3.0 * (__freq - __w3) / (__w0 - __w3)
                            __I1 = 1.0 / 3.0 * (__freq - __w3) / (__w1 - __w3)
                            __I2 = 1.0 / 3.0 * (__freq - __w3) / (__w2 - __w3)
                            __I3 = 1.0 / 3.0 * ((__freq - __w0) / (__w3 - __w0)
                                                + (__freq - __w1) / (__w3 - __w1)
                                                + (__freq - __w2) / (__w3 - __w2))

                            _pdos[:, ind_x] += volume_ith / volume_bz * __g \
                                               * (__I0 * __F0 + __I1 * __F1 + __I2 * __F2 + __I3 * __F3)

                        else:
                            pass

    return _pdos
