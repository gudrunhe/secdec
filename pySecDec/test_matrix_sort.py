"""Unit tests for the algebra module"""

from __future__ import print_function
from .matrix_sort import *
import unittest
import sympy as sp
import numpy as np
from itertools import permutations
from nose.plugins.attrib import attr

class TestMatrixSorting(unittest.TestCase):
    def setUp(self):
        self.matrix_3_by_4 = np.array([[ 2,  2,  1,  1],
                                       [ 2,  0,  2,  1],
                                       [ 2,  2,  0,  0]])

    def check_algorithm(self, sort_function, matrix, solutions):
        for permutation1 in permutations(list(range(3))):
            for permutation2 in permutations(list(range(3))):
                permuted_matrix = matrix
                permuted_matrix[:] = matrix[permutation1,:]
                permuted_matrix[:,1:] = matrix[:,1:][:,permutation2]

                print('permuted')
                print(permuted_matrix)
                sort_function(permuted_matrix)
                print('sorted')
                print(permuted_matrix)
                print('----------------------')
                print()

                found_solution = False
                for solution in solutions:
                    if np.array_equal(permuted_matrix, solution):
                        found_solution = True
                        break
                if not found_solution:
                    raise AssertionError('Did not find an expected solution.')

    #@attr('active')
    def test_iterative_algorithm(self):
        expected_solutions = [
                                 [[2, 0, 0, 2],
                                  [2, 1, 1, 2],
                                  [2, 1, 2, 0]],

                                 [[2, 0, 1, 2],
                                  [2, 2, 0, 0],
                                  [2, 2, 1, 1]]
                             ]
        self.check_algorithm(iterative_sort, self.matrix_3_by_4, expected_solutions)

    #@attr('active')
    def test_light_Pak_algorithm(self):
        expected_solutions = [
                                 [[2, 0, 2, 0],
                                  [2, 1, 2, 1],
                                  [2, 2, 0, 1]],

                                 [[2, 0, 2, 1],
                                  [2, 2, 0, 0],
                                  [2, 2, 1, 1]]
                             ]
        self.check_algorithm(Pak_sort, self.matrix_3_by_4, expected_solutions)

    #@attr('active')
    def test_Pak_algorithm(self):
        expected_solutions = [
                                 [[2, 0, 2, 1],
                                  [2, 2, 0, 0],
                                  [2, 2, 1, 1]]
                             ]
        self.check_algorithm(Pak_sort, self.matrix_3_by_4, expected_solutions)
