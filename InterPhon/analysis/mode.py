import numpy as np
from InterPhon import error
from InterPhon.inout import vasp
from InterPhon.core.pre_check import PreArgument


class Mode(object):
    def __init__(self, process):
        self.process = process
        self.mode_inds = []
        self.k_point = np.empty((3,))
        self.num_images = 30
        # self.freq = np.empty((len(self.process.unit_cell.xyz_true)), dtype=float)
        self.mode = np.empty((len(self.process.unit_cell.xyz_true),
                             len(self.process.unit_cell.xyz_true)))

    def set(self, mode_inds=(0,), k_point=(0.0, 0.0, 0.0)):
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
        _current_position_true = np.transpose(self.process.unit_cell.atom_cart.copy()[self.process.unit_cell.atom_true, :])
        _mass_weight = np.transpose(self.process.unit_cell.mass_true.reshape((-1, 3))) / self.process.unit_cell.mass_true.max()
        for mode_ind in self.mode_inds:
            for ind, x in enumerate(np.linspace(0, 2 * np.pi, self.num_images, endpoint=False), 1):
                self.process.unit_cell.atom_cart[self.process.unit_cell.atom_true, :] = \
                    np.transpose(_current_position_true + np.sin(x)
                                 * np.transpose(self.mode[mode_ind, :].reshape((-1, 3)).real) / np.sqrt(_mass_weight))
                # np.sin(x + np.dot(_q, _current_position_true))

                if ind == 1:
                    self.process.unit_cell.selective = False
                    lines = vasp.write_input_lines(self.process.unit_cell, "unknown system")
                    lines[7] = "Cartesian configuration= %4d" % ind + '\n'
                else:
                    self.process.unit_cell.selective = False
                    line = vasp.write_input_lines(self.process.unit_cell, "unknown system")
                    line[7] = "Cartesian configuration= %4d" % ind + '\n'
                    lines.extend(line[7:])

            with open(out_folder + '/XDATCAR_phonon_mode_{0}_{1}'.format(mode_ind, self.k_point), 'w') as outfile:
                outfile.write("%s" % "".join(lines))

    def plot(self, out_folder='.', unit_cell='POSCAR', code_name='vasp'):
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

    def write_mode_displace(self, out_folder='.', amplitude=5.0):
        # Make file of supercell with displacements along normal mode
        # The displaced supercell along an imaginary mode can be used for structure search
        # test is ongoing
        q = np.dot(self.k_point, self.process.reciprocal_matrix)

        super_cell = self.process.super_cell
        user_arg = self.process.user_arg

        # make an enlargement commensurate with self.k_point
        from fractions import Fraction

        enlargement = np.ones([3, ], dtype=int)
        for ind, value in enumerate(user_arg):
            if value:
                enlargement[ind] = Fraction(self.k_point[ind]).limit_denominator(100).denominator
        user_arg.enlargement = enlargement

        super_cell.initialization()
        super_cell.set_super_cell(self.process.unit_cell, user_arg)
        super_cell.set_super_ind_true(self.process.unit_cell, user_arg)
        super_cell.set_mass_true()

        _current_position_true = np.transpose(super_cell.atom_cart.copy()[super_cell.atom_true, :])
        _mass_weight = np.transpose(super_cell.mass_true.reshape((-1, 3))) / super_cell.mass_true.max()
        for mode_ind in self.mode_inds:
            super_cell.atom_cart[super_cell.atom_true, :] = \
                np.transpose(_current_position_true + amplitude * np.sin(np.dot(q, _current_position_true))
                             * self.mode[mode_ind, :].reshape((-1, 3)).real / np.sqrt(_mass_weight))

            _lines = vasp.write_input_lines(self, 'Commensurate supercell with displacements along normal mode {0} at k-point {1}'.format(mode_ind, self.k_point))
            with open(out_folder + '/MPOSCAR-{0}'.format(mode_ind), 'w') as outfile:
                outfile.write("%s" % "".join(_lines))
