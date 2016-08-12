"""
Matrix Sort
-----------

Algorithms to sort a matrix when column and
row permutations are allowed.

"""

from .misc import argsort_2D_array, argsort_ND_array
import numpy as np

def iterative_sort(matrix):
    '''
    Inplace modify the `matrix` to some ordering,
    when permutations of rows and columns (excluding
    the first) are allowed.

    .. note::
        This function may result in different
        orderings depending on the initial ordering.

    .. seealso::
        :func:`.Pak_sort`

    :param matrix:
        2D array-like;
        The matrix to be canonicalized.

    '''
    # keep sorting the 2D matrix (in place) along both dimensions
    # until it no longer changes

    # initilize sort keys such that we enter the while loop at least once
    sort_key_axis_0 = None
    sort_key_axis_1 = None

    while not ( np.array_equal(sort_key_axis_0, np.arange(matrix.shape[0])) and np.array_equal(sort_key_axis_1, np.arange(matrix.shape[1]-1)) ):
            # sort along axis 0
            sort_key_axis_0 = argsort_2D_array(matrix)
            matrix[:] = matrix[sort_key_axis_0]

            # sort along axis 1 excluding the coefficients (first column)
            sort_key_axis_1 = argsort_2D_array(matrix.T[1:])
            matrix[:,1:] = matrix[:,1:][:,sort_key_axis_1]

def Pak_sort(matrix):
    '''
    Inplace modify the `matrix` to some ordering,
    when permutations of rows and columns (excluding
    the first) are allowed. The implementation of
    this function is describedin chapter 2 of
    [Pak11]_.

    .. note::
        This function may result in different
        orderings depending on the initial ordering.

    .. seealso::
        :func:`.iterative_sort`

    :param matrix:
        2D array-like;
        The matrix to be canonicalized.

    '''
    for i in range(1,matrix.shape[1]):
        # sort all permutations of columns `i` and `j` (where `i`<=`j`) --> pick largest
        permutaions = []
        for j in range(i,matrix.shape[1]):
            permuted_matrix = matrix.copy()

            # permute integration variables `i` and `j`
            permuted_matrix[:,i] = matrix[:,j]
            permuted_matrix[:,j] = matrix[:,i]

            # sort by rows
            permuted_matrix[:] = permuted_matrix[argsort_2D_array(matrix)]

            # transpose since we need column-wise ordering in the next step
            permutaions.append(permuted_matrix.T)

        # find the largest `i`th column in `permutaions` and keep that
        index_of_largest = argsort_ND_array(permutaions)[-1]
        matrix[:] = permutaions[index_of_largest].T # transpose back
