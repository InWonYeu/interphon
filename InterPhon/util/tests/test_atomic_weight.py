import numpy as np
import unittest

from InterPhon.util import get_atomic_weight


class TestAtomicWeight(unittest.TestCase):
    def test_get_atomic_weight(self):
        mass_H = get_atomic_weight('H')

        self.assertAlmostEqual(mass_H, 1.0, delta=1e-02)
        with self.assertRaises(NameError):
            get_atomic_weight('Ha')


if __name__ == "__main__":
    unittest.main()
