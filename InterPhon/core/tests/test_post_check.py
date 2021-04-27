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
                                   atom_cart=np.array([[0, 0, 0], [0.5, 0, 0], [0.25, 0.25, 0.25], [0.75, 0.25, 0.25]]),
                                   atom_true=[0, 1],
                                   xyz_true=[0, 0, 0, 1, 1, 1],
                                   mass_true=None)

    def test_displacement(self):
        user_args_ = PostArgument()
        user_args_.displacement = self.displacement

        self.assertIsInstance(user_args_.displacement, float)
        self.assertEqual(user_args_.displacement, 0.02)

    def test_enlargement(self):
        user_args_ = PostArgument()
        user_args_.enlargement = self.enlargement

        self.assertIsInstance(user_args_.enlargement, np.ndarray)
        self.assertListEqual(list(user_args_.enlargement), [2, 2, 1])
        with self.assertRaises(error.Insufficient_ENLARGE_Error):
            user_args_.enlargement = [2, 2]

    def test_periodicity(self):
        user_args_ = PostArgument()
        user_args_.enlargement = self.enlargement
        user_args_.periodicity = self.periodicity

        self.assertIsInstance(user_args_.periodicity, np.ndarray)
        self.assertListEqual(list(user_args_.periodicity), [1, 1, 0])
        with self.assertRaises(error.Insufficient_PBC_Error):
            user_args_.periodicity = "True True"

        with self.assertRaises(error.Mismatch_ENLARGE_and_PBC_Error):
            user_args_.periodicity = "True False True"

    def test_set_user_argument(self):
        user_args_ = PostArgument()
        user_args_.set_user_argument(self.user_args)

        self.assertEqual(user_args_.displacement, 0.02)
        self.assertListEqual(list(user_args_.enlargement), [2, 2, 1])
        self.assertListEqual(list(user_args_.periodicity), [1, 1, 0])

        user_args_.initialization()
        self.assertIsNone(user_args_.displacement)
        self.assertIsNone(user_args_.enlargement)
        self.assertIsNone(user_args_.periodicity)

    def test_check_user_argument(self):
        user_args_ = PostArgument()
        user_args_.set_user_argument(self.user_args)
        user_args_.enlargement = [2, 2, 2]

        with self.assertRaises(error.Mismatch_ENLARGE_and_PBC_Error):
            user_args_.check_user_argument()

    def test_check_match_argument(self):
        user_args_ = PostArgument()
        user_args_.set_user_argument(self.user_args)

        with self.assertRaises(error.Mismatch_ENLARGE_post_Error):
            user_args_.check_match_argument(self.unit_cell, self.super_cell)


if __name__ == "__main__":
    unittest.main()
