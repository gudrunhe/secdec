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

                                 [[2, 0, 0, 2],
                                  [2, 1, 1, 2],
                                  [2, 1, 2, 0]]
                             ]
        self.check_algorithm(light_Pak_sort, self.matrix_3_by_4, expected_solutions)

    #@attr('active')
    def test_Pak_algorithm(self):
        expected_solutions = [
                                 [[2, 0, 0, 2],
                                  [2, 1, 1, 2],
                                  [2, 1, 2, 0]]
                             ]
        self.check_algorithm(Pak_sort, self.matrix_3_by_4, expected_solutions)

    def test_Pak_indices(self):

        expolist_initial = np.array([ [4, 3],
                                      [2, 1] ])

        # Canonical with start_column = 0
        target_expolist_canonical0 = np.array([ [1, 2],
                                                [3, 4] ])
        # Canonical with start_column = 1
        target_expolist_canonical1 = np.array([ [2, 1],
                                                [4, 3] ])

        expolist_canonical0 = expolist_initial.copy()
        Pak_sort(expolist_canonical0,[0,1])
        np.testing.assert_array_equal(expolist_canonical0,target_expolist_canonical0)

        expolist_canonical1 = expolist_initial.copy()
        Pak_sort(expolist_canonical1)
        np.testing.assert_array_equal(expolist_canonical1,target_expolist_canonical1)

    def test_Pak_indices_groups(self):

        target = np.array([ [0,1 , 0,2],
                            [1,0 , 1,0] ])

        matrices1 = [ [ [1,0],[0,1] ],
                      [ [0,1],[1,0] ] ]

        matrices2 =[ [ [1,0],[0,2] ],
                     [ [0,1],[2,0] ],
                     [ [2,0],[0,1] ],
                     [ [0,2],[1,0] ] ]

        for mat1 in matrices1:
            for mat2 in matrices2:
                permutation = np.hstack( (mat1,mat2) )
                Pak_sort(permutation, [0, 1], [2, 3])
                np.testing.assert_array_equal(permutation, target)
