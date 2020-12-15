import numpy as np
from InterPhon import error
from InterPhon.inout import vasp


class Mode(object):
    def __init__(self, process):
        self.process = process
        self.mode_inds = []
        self.k_point = np.empty((3,))
        self.num_images = 30
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
            raise error.Not_Specified_Kpath_Error

        self.mode = self.process.v_q[_ind_k, :, :]

    def write(self, out_folder='.'):
        _current_position_true = self.process.unit_cell.atom_cart.copy()[self.process.unit_cell.atom_true]
        _mass_weight = self.process.unit_cell.mass_true.reshape((-1, 3)) / self.process.unit_cell.mass_true.max()
        for mode_ind in self.mode_inds:
            for ind, x in enumerate(np.linspace(0, 2 * np.pi, self.num_images, endpoint=False), 1):
                self.process.unit_cell.atom_cart[self.process.unit_cell.atom_true] = _current_position_true \
                    + np.sin(x) * self.mode[mode_ind, :].reshape((-1, 3)).real / np.sqrt(_mass_weight)

                if ind == 1:
                    lines = vasp.write_input_lines(self.process.unit_cell, "unknown system")
                    lines[7] = "Cartesian configuration= %4d" % ind + '\n'
                else:
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
        _current_position_true = atom.positions.copy()[self.process.unit_cell.atom_true]
        _mass_weight = self.process.unit_cell.mass_true.reshape((-1, 3)) / self.process.unit_cell.mass_true.max()
        for mode_ind in self.mode_inds:
            traj = Trajectory(out_folder + "/Trajectory_{0}.traj".format(mode_ind), 'w')
            for ind, x in enumerate(np.linspace(0, 2 * np.pi, self.num_images, endpoint=False), 1):
                atom.positions[self.process.unit_cell.atom_true] = _current_position_true \
                    + np.sin(x) * self.mode[mode_ind, :].reshape((-1, 3)).real / np.sqrt(_mass_weight)
                traj.write(atom)
            traj.close()
            atoms = iread(out_folder + "/Trajectory_{0}.traj".format(mode_ind))
            view(atoms)
