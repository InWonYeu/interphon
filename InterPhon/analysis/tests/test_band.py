import numpy as np
import unittest

from InterPhon.analysis import band


class TestBand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x = tf.random.normal(shape=(1, 3, 4))
        cls.index = tf.convert_to_tensor(value=[0, 0, 0, 1, 1, 2])
        cls.n = tf.convert_to_tensor(value=[3, 2, 1])

    def test_set(self):
        repeat_result = _repeat(self.x, self.n, axis=1).numpy()
        self.assertListEqual(list(repeat_result.shape), [1, 6, 4])

    def test_write(self):
        repeat_result = _repeat(self.x, self.n, axis=1).numpy()
        self.assertEqual(repeat_result[0, 0, 0], repeat_result[0, 1, 0])

    def test_plot(self):
        repeat_result = _repeat(self.x, self.n, axis=1).numpy()
        self.assertNotEqual(repeat_result[0, 0, 0], repeat_result[0, 3, 0])

    def test_plot_with_dos(self):
        repeat_result = _repeat(self.x, self.n, axis=1).numpy()
        self.assertEqual(repeat_result[0, 3, 0], repeat_result[0, 4, 0])


if __name__ == "__main__":
    unittest.main()
