import numpy as np
import unittest

from InterPhon.core import PostArgument, UnitCell, SuperCell
from InterPhon import error


class TestPreArgument(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.displacement = 0.02
        cls.enlargement = (2, 2, 1)
        cls.periodicity = "True True False"
        cls.user_args = {'displacement': 0.02, 'enlargement': "2 2 1", 'periodicity': "1 1 0"}

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

    def test_check_match_argument(self):
        user_args_ = PostArgument()
        user_args_.set_user_argument(self.user_args)

        with self.assertRaises(error.Mismatch_ENLARGE_post_Error):
            user_args_.check_match_argument(self.unit_cell, self.super_cell)


if __name__ == "__main__":
    unittest.main()
