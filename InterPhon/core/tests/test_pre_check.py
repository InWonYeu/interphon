import numpy as np
import unittest

from InterPhon.core import PreArgument
from InterPhon import error


class TestPreArgument(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.displacement = 0.02
        cls.enlargement = (2, 2, 1)
        cls.periodicity = "True True False"
        cls.user_args = {'displacement': 0.02, 'enlargement': "2 2 1", 'periodicity': "1 1 0"}

    def test_displacement(self):
        user_args_ = PreArgument()
        user_args_.displacement = self.displacement

        self.assertIsInstance(user_args_.displacement, float)
        self.assertEqual(user_args_.displacement, 0.02)

    def test_enlargement(self):
        user_args_ = PreArgument()
        user_args_.enlargement = self.enlargement

        self.assertIsInstance(user_args_.enlargement, np.ndarray)
        self.assertListEqual(list(user_args_.enlargement), [2, 2, 1])
        with self.assertRaises(error.Insufficient_ENLARGE_Error):
            user_args_.enlargement = [2, 2]

    def test_periodicity(self):
        user_args_ = PreArgument()
        user_args_.enlargement = self.enlargement
        user_args_.periodicity = self.periodicity

        self.assertIsInstance(user_args_.periodicity, np.ndarray)
        self.assertListEqual(list(user_args_.periodicity), [1, 1, 0])
        with self.assertRaises(error.Insufficient_PBC_Error):
            user_args_.periodicity = "True True"

        with self.assertRaises(error.Mismatch_ENLARGE_and_PBC_Error):
            user_args_.periodicity = "True False True"

    def test_set_user_argument(self):
        user_args_ = PreArgument()
        user_args_.set_user_argument(self.user_args)

        self.assertEqual(user_args_.displacement, 0.02)
        self.assertListEqual(list(user_args_.enlargement), [2, 2, 1])
        self.assertListEqual(list(user_args_.periodicity), [1, 1, 0])

        user_args_.initialization()
        self.assertIsNone(user_args_.displacement)
        self.assertIsNone(user_args_.enlargement)
        self.assertIsNone(user_args_.periodicity)

    def test_check_user_argument(self):
        user_args_ = PreArgument()
        user_args_.set_user_argument(self.user_args)
        user_args_.enlargement = [2, 2, 2]

        with self.assertRaises(error.Mismatch_ENLARGE_and_PBC_Error):
            user_args_.check_user_argument()


if __name__ == "__main__":
    unittest.main()
