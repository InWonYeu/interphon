import numpy as np
from InterPhon import error
from InterPhon.inout import vasp, aims, espresso
from InterPhon.util import MatrixLike


def to_int_numpy(data):
    try:
        if isinstance(data, int):
            data = [data, ]
        elif isinstance(data, str):
            data = [int(val) for val in data.strip().split()]
        else:
            data = [int(val) for val in data]
    except ValueError:
        raise ValueError("The items of '{0}' cannot be converted to int".format(data))
    return np.array(data)


class Mode(object):
    """
    Mode class to analyze the phonon dispersion.
    Its instance variables contain an instance of :class:`core.PostProcess` class storing the information on eigen-frequency and eigen-mode.
    Based on the given information, phonon modes are visualized at a specified k-point.
    The vibrational motions are written to external data and movie files using the :class:`analysis.Mode.write`
    and :class:`analysis.Mode.plot` methods, respectively.

    :param process: Instance of PostProcess class
    :type process: :class:`core.PostProcess`
    """
    def __init__(self, process):
        """
        Constructor of Mode class.
        """
        self.process = process
        self.__mode_inds = []
        self.k_point = np.empty((3,))
        self.num_images = 30
        # self.freq = np.empty((len(self.process.unit_cell.xyz_true)), dtype=float)
        self.mode = np.empty((len(self.process.unit_cell.xyz_true),
                             len(self.process.unit_cell.xyz_true)))

    @property
    def mode_inds(self):
        return self.__mode_inds

    @mode_inds.setter
    def mode_inds(self, _mode_inds):
        _mode_inds = to_int_numpy(_mode_inds)
        self.__mode_inds = _mode_inds

    def set(self, mode_inds: MatrixLike = (0,),
            k_point: MatrixLike = (0.0, 0.0, 0.0)):
        """
        Set the k-point, mode index, and corresponding phonon mode.

        :param mode_inds: The index of phonon mode labeled by (k-point, index)
        :type mode_inds: MatrixLike or int
        :param k_point: The k-point of phonon mode labeled by (k-point, index)
        :type k_point: MatrixLike
        """
        self.mode_inds = mode_inds
        self.k_point = np.array(k_point)

        _ind_k = None
        for ind, kpt in enumerate(self.process.k_points):
            if np.allclose(self.k_point, kpt):
                _ind_k = ind
        if _ind_k is None:
            raise error.Not_Specified_Kpath_Error(self.k_point)

        # self.freq = self.process.w_q[_ind_k, :]
        self.mode = self.process.v_q[_ind_k, :, :]

    def write(self, out_folder='.'):
        """
        Write phonon Mode in the file name of **XDATCAR_phonon_mode**.

        :param out_folder: Folder path for **XDATCAR_phonon_mode** to be stored, defaults to .
        :type out_folder: str
        """
        from copy import deepcopy

        unit_cell = deepcopy(self.process.unit_cell)
        unit_cell.selective = False

        _current_position_true = np.transpose(unit_cell.atom_cart.copy()[unit_cell.atom_true, :])
        _mass_weight = np.transpose(unit_cell.mass_true.reshape((-1, 3))) / unit_cell.mass_true.max()
        for mode_ind in self.mode_inds:
            for ind, x in enumerate(np.linspace(0, 2 * np.pi, self.num_images, endpoint=False), 1):
                unit_cell.atom_cart[unit_cell.atom_true, :] = \
                    np.transpose(_current_position_true + np.sin(x)
                                 * np.transpose(self.mode[mode_ind, :].reshape((-1, 3)).real) / np.sqrt(_mass_weight))
                # np.sin(x + np.dot(_q, _current_position_true))

                if ind == 1:
                    lines = vasp.write_input_lines(unit_cell, "unknown system")
                    lines[7] = "Cartesian configuration= %4d" % ind + '\n'
                else:
                    line = vasp.write_input_lines(unit_cell, "unknown system")
                    line[7] = "Cartesian configuration= %4d" % ind + '\n'
                    lines.extend(line[7:])

            with open(out_folder + '/XDATCAR_phonon_mode_{0}_{1}'.format(mode_ind, self.k_point), 'w') as outfile:
                outfile.write("%s" % "".join(lines))

    def plot(self, out_folder='.',
             unit_cell='POSCAR',
             code_name='vasp'):
        """
        Visualize phonon Mode using modules of the `Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/index.html>`_.

        :param out_folder: Folder path for **Trajectory.traj** to be stored, defaults to .
        :type out_folder: str
        :param unit_cell: Path of unit cell input file, defaults to POSCAR
        :type unit_cell: str
        :param code_name: Specification of the file-format by a DFT program, defaults to vasp
        :type code_name: str
        """
        try:
            from ase.io.trajectory import Trajectory
            from ase.io import read, iread
            from ase.visualize import view
        except ImportError:
            raise ImportError("\nThe parent directory of ase package must be included in 'sys.path'")

        if code_name == 'espresso':
            code_name = 'espresso-in'  # aims, espresso-in, vasp

        atom = read(unit_cell, format=code_name)
        _current_position_true = np.transpose(atom.positions.copy()[self.process.unit_cell.atom_true, :])
        _mass_weight = np.transpose(self.process.unit_cell.mass_true.reshape((-1, 3))) / self.process.unit_cell.mass_true.max()
        for mode_ind in self.mode_inds:
            traj = Trajectory(out_folder + "/Trajectory_{0}.traj".format(mode_ind), 'w')
            for _, x in enumerate(np.linspace(0, 2 * np.pi, self.num_images, endpoint=False), 1):
                atom.positions[self.process.unit_cell.atom_true, :] = \
                    np.transpose(_current_position_true + np.sin(x)
                                 * np.transpose(self.mode[mode_ind, :].reshape((-1, 3)).real) / np.sqrt(_mass_weight))
                # np.sin(x + np.dot(_q, _current_position_true))

                traj.write(atom)

            traj.close()
            atoms = iread(out_folder + "/Trajectory_{0}.traj".format(mode_ind))
            view(atoms)

    def write_mode_displace(self, out_folder='.',
                            amplitude: float =1.0,
                            code_name: str = 'vasp'):
        """
        Write a supercell file which is commensurate with k-point and displaced along phonon Mode in the file name of **MPOSCAR**.

        :param out_folder: Folder path for **MPOSCAR** to be stored, defaults to .
        :type out_folder: str
        :param amplitude: Amplitude of vibrational motions, defaults to 1.0
        :type amplitude: float
        :param code_name: Specification of the file-format by a DFT program, defaults to vasp
        :type code_name: str
        """
        # Make a supercell file with displacement along normal mode
        # The displaced supercell along an imaginary mode can be used for structure search
        # test is ongoing
        from copy import deepcopy
        from fractions import Fraction

        q = np.dot(self.k_point, self.process.reciprocal_matrix)
        super_cell = deepcopy(self.process.super_cell)
        user_arg = deepcopy(self.process.user_arg)

        # make an enlargement commensurate with self.k_point
        enlargement = np.ones([3, ], dtype=int)
        for ind, value in enumerate(user_arg.periodicity):
            if value:
                enlargement[ind] = Fraction(self.k_point[ind]).limit_denominator(100).denominator
        user_arg.enlargement = enlargement

        _enlarge = 1
        for ind, value in enumerate(user_arg.periodicity):
            if value:
                _enlarge = _enlarge * user_arg.enlargement[ind]

        super_cell.initialization()
        super_cell.set_super_cell(self.process.unit_cell, user_arg)
        super_cell.set_super_ind_true(self.process.unit_cell, user_arg)
        super_cell.set_mass_true()
        super_cell.selective = True

        _current_position_true = np.transpose(super_cell.atom_cart.copy()[super_cell.atom_true, :])
        _mass_weight = super_cell.mass_true.reshape((-1, 3)) / super_cell.mass_true.max()
        _phase_factor = np.cos(np.dot(q, _current_position_true).reshape((-1, 1)))
        for mode_ind in self.mode_inds:
            mode_super_cell = []
            for mode in self.mode[mode_ind, :].reshape((-1, 3)).real:
                for i in range(_enlarge):
                    mode_super_cell.append(mode)
            mode_super_cell = np.asfarray(mode_super_cell)

            super_cell.atom_cart[super_cell.atom_true, :] = np.transpose(_current_position_true) + amplitude \
                                                            * _phase_factor * mode_super_cell / np.sqrt(_mass_weight)

            comment = 'Commensurate supercell with displacements along normal mode {0} at k-point {1}'.format(mode_ind, self.k_point)
            if code_name == 'vasp':
                _lines = vasp.write_input_lines(super_cell, comment)
                with open(out_folder + '/MPOSCAR-{0}'.format(mode_ind), 'w') as outfile:
                    outfile.write("%s" % "".join(_lines))

            elif code_name == 'espresso':
                _lines = espresso.write_input_lines(super_cell, comment)
                with open(out_folder + '/MPOSCAR-{0}'.format(mode_ind), 'w') as outfile:
                    outfile.write("%s" % "".join(_lines))

            elif code_name == 'aims':
                _lines = aims.write_input_lines(super_cell, comment)
                with open(out_folder + '/MPOSCAR-{0}'.format(mode_ind), 'w') as outfile:
                    outfile.write("%s" % "".join(_lines))
