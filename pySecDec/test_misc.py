"""Unit tests for the miscellaneous routines"""

from .misc import *
import numpy as np
import unittest

class TestSort(unittest.TestCase):
    def test_sort_2D_array(self):
        in_array = [[1,2,3],
                    [2,3,4],
                    [1,2,3]]
        target_sorted_array = [[1,2,3],
                               [1,2,3],
                               [2,3,4]]
        calculated_sort_indices = argsort_2D_array(in_array)

        in_array = np.asarray(in_array)
        np.testing.assert_array_equal(in_array[calculated_sort_indices], target_sorted_array)

    def test_error_message(self):
        self.assertRaisesRegexp(AssertionError, 'array.*two dimensional', argsort_2D_array, [1,2,3,4])

class TestPowerset(unittest.TestCase):
    def test_powerset_range(self):
        self.assertEqual(list(powerset(range(0))), [()])
        self.assertEqual(list(powerset(range(1))), [(),(0,)])
        self.assertEqual(list(powerset(range(2))), [(),(0,),(1,),(0,1)])
        self.assertEqual(list(powerset(range(3))), [(),(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)])
        self.assertEqual(list(powerset(range(4))), [(),(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)])

    def test_powerset_without_empty(self):
        self.assertEqual(list(powerset(range(0),exclude_empty=True)), [])
        self.assertEqual(list(powerset(range(1),exclude_empty=True)), [(0,)])
        self.assertEqual(list(powerset(range(2),exclude_empty=True)), [(0,),(1,),(0,1)])
        self.assertEqual(list(powerset(range(3),exclude_empty=True)), [(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)])
        self.assertEqual(list(powerset(range(4),exclude_empty=True)), [(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)])

class TestDet(unittest.TestCase):
    def test_calculation(self):
        M1 = [[1,1],
              [2,2]]
        target_det_M1 = 0
        self.assertEqual(det(M1), target_det_M1)

        M2 = [[1,2,3],
              [0,2,2],
              [0,0,2]]
        target_det_M2 = 4
        self.assertEqual(det(M2), target_det_M2)

    def test_error_messages(self):
        # wrong shape
        M1 = [1,2,3]
        M2 = [[[1,2],[3,4]],[[5,6],[7,8]],[[9,10],[11,12]]]
        for M in (M1,M2):
            self.assertRaisesRegexp(AssertionError, 'M.*two dimensional', det, M)

        # not square
        M3 = [[1,2,3],
              [0,2,2]]
        self.assertRaisesRegexp(AssertionError, 'M.*square', det, M3)

    def test_1x1_special_case(self):
        M = [[1]]
        self.assertEqual(det(M), 1)
        self.assertTrue(np.issubdtype(type(det(M)), np.int))

class TestAdjugate(unittest.TestCase):
    def test_calculation(self):
        M = [[1,2,3],
             [4,5,6],
             [7,8,9]]
        target_adj_M = [[-3,   6, -3],
                        [ 6, -12,  6],
                        [-3,   6, -3]]
        np.testing.assert_array_equal(adjugate(M), target_adj_M)

    def test_error_messages(self):
        # wrong shape
        M1 = [1,2,3]
        M2 = [[[1,2],[3,4]],[[5,6],[7,8]],[[9,10],[11,12]]]
        for M in (M1,M2):
            self.assertRaisesRegexp(AssertionError, 'M.*two dimensional', adjugate, M)

        # not square
        M3 = [[1,2,3],
              [0,2,2]]
        self.assertRaisesRegexp(AssertionError, 'M.*square', adjugate, M3)

    def test_1x1_special_case(self):
        M1 = [[1 ]]
        M2 = [[5 ]]
        M3 = [[2.]]
        for M in (M1,M2,M3):
            adj_M = adjugate(M)
            self.assertEqual(adj_M, 1)
            self.assertTrue(isinstance(adj_M, np.ndarray))

        # should keep data type of input
        self.assertTrue(np.issubdtype(adjugate(M1).dtype, np.int))
        self.assertTrue(np.issubdtype(adjugate(M2).dtype, np.int))
        self.assertTrue(np.issubdtype(adjugate(M3).dtype, np.float))
