import numpy as np
import unittest

from InterPhon.core import PreArgument, UnitCell, SuperCell
from InterPhon import error


class TestPreArgument(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.user_args_ = {'displacement': 0.02, 'enlargement': "2 1 1", 'periodicity': "1 0 0"}
        cls.user_args = PreArgument()
        cls.user_args.set_user_argument(cls.user_args_)

        cls.unit_cell = UnitCell(lattice_matrix=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                                 atom_type=['H', 'He'],
                                 num_atom=np.array([1, 1]),
                                 coordinate='cartesian',
                                 atom_cart=np.array([[0, 0, 0], [0.5, 0.5, 0.5]]),
                                 atom_true=[0],
                                 xyz_true=[0, 0, 0],
                                 mass_true=None)

        cls.super_cell = SuperCell(lattice_matrix=np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
                                   atom_type=['H', 'H', 'He', 'He'],
                                   num_atom=np.array([2, 2]),
                                   coordinate='cartesian',
                                   atom_cart=np.array([[0, 0, 0], [1.0, 0, 0], [0.5, 0.5, 0.5], [1.5, 0.5, 0.5]]),
                                   atom_true=[0, 1],
                                   xyz_true=[0, 0, 0, 1, 1, 1],
                                   mass_true=None)
        cls.super_cell.set_mass_true()

    def test_set_super_cell(self):
        super_cell_ = SuperCell()
        super_cell_.set_super_cell(unit_cell=self.unit_cell, user_arg=self.user_args)

        self.assertListEqual([list(v) for v in super_cell_.lattice_matrix],
                             [list(v) for v in self.super_cell.lattice_matrix])
        self.assertListEqual(super_cell_.atom_type, self.super_cell.atom_type)
        self.assertListEqual(list(super_cell_.num_atom), list(self.super_cell.num_atom))
        self.assertEqual(super_cell_.coordinate, self.super_cell.coordinate)
        self.assertListEqual([list(v) for v in super_cell_.atom_cart],
                             [list(v) for v in self.super_cell.atom_cart])

    def test_set_super_ind_true(self):
        super_cell_ = SuperCell()
        super_cell_.set_super_cell(unit_cell=self.unit_cell, user_arg=self.user_args)
        super_cell_.set_super_ind_true(unit_cell=self.unit_cell, user_arg=self.user_args)

        self.assertListEqual(super_cell_.atom_true, self.super_cell.atom_true)
        self.assertListEqual(super_cell_.xyz_true, self.super_cell.xyz_true)

        super_cell_.set_mass_true()

        self.assertIsNotNone(super_cell_.mass_true)
        self.assertListEqual(list(super_cell_.mass_true), list(self.super_cell.mass_true))


if __name__ == "__main__":
    unittest.main()
