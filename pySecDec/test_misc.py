"""Unit tests for the miscellaneous routines"""

from .misc import *
import numpy as np
import unittest
from nose.plugins.attrib import attr

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

    #@attr('active')
    def test_strided_powerset(self):
        self.assertEqual(list(powerset(range(4),stride=2)), [(),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2,3)])
        self.assertEqual(list(powerset(range(4),stride=3)), [(),(0,1,2),(0,1,3),(0,2,3),(1,2,3)])
        self.assertEqual(list(powerset(range(4),stride=2,exclude_empty=True)), [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2,3)])

class TestMissing(unittest.TestCase):
    #@attr('active')
    def test_missing(self):
        self.assertEqual(missing([1,2,3], [1,2]), [3])
        self.assertEqual(missing(['a','b','c','d','e','f'], ['a','e','d']), ['b','c','f'])

    #@attr('active')
    def test_in_combination_with_powerset(self):
        full_set = list(range(4))
        powerset_4_generator = powerset(full_set)
        target_powerset_4 = [(),(0,),(1,),(2,),(3,),(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,1,2,3)]

        target_missing_in_powerset_4 = list(target_powerset_4)
        target_missing_in_powerset_4.reverse()

        for i, powerset_item in enumerate(powerset_4_generator):
            print(i)
            self.assertEqual(powerset_item, target_powerset_4[i])
            self.assertEqual(missing(full_set, powerset_item), list(target_missing_in_powerset_4[i]))

class TestAllPairs(unittest.TestCase):
    #@attr('active')
    def test_error_input_not_even(self):
        self.assertRaisesRegexp(AssertionError, 'even', all_pairs, [1,2,3])

    #@attr('active')
    def test_partitioning(self):
        # list input
        lst = [1,2,3,4]
        self.assertEqual(list(all_pairs(lst)), [[(1,2),(3,4)], [(1,3),(2,4)], [(1,4),(2,3)]])

        # iterator input
        generator = (i for i in lst)
        self.assertEqual(list(all_pairs(generator)), [[(1,2),(3,4)], [(1,3),(2,4)], [(1,4),(2,3)]])

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

class TestCachedProperty(unittest.TestCase):
    #@attr('active')
    def test_usability(self):
        class C(object):
            'Sum up the numbers from one to `N`.'
            def __init__(self, N):
                self.N = N
            @cached_property
            def sum(self):
                result = 0
                for i in range(1, self.N + 1):
                    result += i
                return result
        to_ten = C(10)
        self.assertEqual(to_ten.sum, 1+2+3+4+5+6+7+8+9+10)

    #@attr('active')
    def test_called_only_once(self):
        class C(object):
            def __init__(self):
                self.counter = 0
            @cached_property
            def prop(self):
                self.counter += 1
                return 5 + 5
        instance = C()

        self.assertEqual(instance.counter, 0)
        for i in range(10):
            # the method should be called only once
            self.assertEqual(instance.prop, 10)
            self.assertEqual(instance.counter, 1)
