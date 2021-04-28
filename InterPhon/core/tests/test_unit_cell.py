import numpy as np
import unittest

from InterPhon.core import UnitCell
from InterPhon import error


class TestUnitCell(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lattice_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        cls.atom_type = ['H', 'He']
        cls.num_atom = [1, 1]
        cls.coordinate = 'cartesian'
        cls.atom_cart = [[0, 0, 0], [0.5, 0.5, 0.5]]
        cls.atom_true = [0]
        cls.xyz_true = [0, 0, 0]

    def test_lattice_matrix(self):
        unit_cell_ = UnitCell()
        unit_cell_.lattice_matrix = self.lattice_matrix

        self.assertIsInstance(unit_cell_.lattice_matrix, np.ndarray)
        self.assertListEqual([list(v) for v in unit_cell_.lattice_matrix], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        with self.assertRaises(ValueError):
            unit_cell_.lattice_matrix = [[1, 0, 0], [0, 1, 0], [0, 0, ]]

    def test_atom_type(self):
        unit_cell_ = UnitCell()
        unit_cell_.atom_cart = self.atom_cart
        unit_cell_.atom_type = self.atom_type

        self.assertIsInstance(unit_cell_.atom_type, list)
        for _atom_type in unit_cell_.atom_type:
            self.assertIsInstance(_atom_type, str)

        self.assertListEqual(unit_cell_.atom_type, ['H', 'He'])
        with self.assertRaises(ValueError):
            unit_cell_.atom_type = ['H', ]

    def test_num_atom(self):
        unit_cell_ = UnitCell()
        unit_cell_.atom_cart = self.atom_cart
        unit_cell_.num_atom = self.num_atom

        self.assertIsInstance(unit_cell_.num_atom, np.ndarray)
        self.assertListEqual(list(unit_cell_.num_atom), [1, 1])
        with self.assertRaises(ValueError):
            unit_cell_.num_atom = [1, ]

    def test_coordinate(self):
        unit_cell_ = UnitCell()
        unit_cell_.coordinate = self.coordinate.title()

        self.assertIsInstance(unit_cell_.coordinate, str)
        self.assertEqual(unit_cell_.coordinate, 'cartesian')
        with self.assertRaises(ValueError):
            unit_cell_.coordinate = 'Wrong coordinate'

    def test_atom_cart(self):
        unit_cell_ = UnitCell()
        unit_cell_.atom_cart = self.atom_cart
        unit_cell_.atom_type = self.atom_type

        self.assertIsInstance(unit_cell_.atom_cart, np.ndarray)
        self.assertListEqual([list(v) for v in unit_cell_.atom_cart], [[0, 0, 0], [0.5, 0.5, 0.5]])
        with self.assertRaises(ValueError):
            unit_cell_.atom_cart = [[0, 0, 0], ]

    def test_atom_true(self):
        unit_cell_ = UnitCell()
        unit_cell_.atom_type = self.atom_type
        unit_cell_.atom_true = self.atom_true

        self.assertIsInstance(unit_cell_.atom_true, list)
        for _atom_true in unit_cell_.atom_true:
            self.assertIsInstance(_atom_true, int)

        self.assertListEqual(unit_cell_.atom_true, [0])
        with self.assertRaises(ValueError):
            unit_cell_.atom_true = [2]

    def test_xyz_true(self):
        unit_cell_ = UnitCell()
        unit_cell_.atom_true = self.atom_true
        unit_cell_.xyz_true = self.xyz_true

        self.assertIsInstance(unit_cell_.xyz_true, list)
        for _xyz_true in unit_cell_.xyz_true:
            self.assertIsInstance(_xyz_true, int)

        self.assertListEqual(unit_cell_.xyz_true, [0, 0, 0])
        with self.assertRaises(ValueError):
            unit_cell_.xyz_true = [2, 0, 0]

    def test_initialization(self):
        unit_cell_ = UnitCell()
        unit_cell_.lattice_matrix = self.lattice_matrix
        unit_cell_.atom_type = self.atom_type
        unit_cell_.num_atom = self.num_atom
        unit_cell_.coordinate = self.coordinate
        unit_cell_.atom_cart = self.atom_cart
        unit_cell_.atom_true = self.atom_true
        unit_cell_.xyz_true = self.xyz_true

        unit_cell_.set_mass_true()
        self.assertIsInstance(unit_cell_.mass_true, np.ndarray)
        self.assertListEqual(list(unit_cell_.mass_true.shape), [3 * len(self.atom_true), ])

        unit_cell_.initialization()
        for key, value in unit_cell_.__dict__.items():
            self.assertIsNone(value)
