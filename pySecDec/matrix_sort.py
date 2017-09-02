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
        :func:`.Pak_sort`, :func:`.light_Pak_sort`

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
    this function is described in chapter 2 of
    [Pak11]_.

    .. note::
        This function may result in different
        orderings depending on the initial ordering.

    .. seealso::
        :func:`.iterative_sort`, :func:`.light_Pak_sort`

    :param matrix:
        2D array-like;
        The matrix to be canonicalized.

    '''
    options = [matrix]
    for i in range(1, matrix.shape[1]):
        permutations = []
        for m in options:
            for j in range(i, m.shape[1]):
                permuted_matrix = m.copy()

                # permute integration variables `column_to_swap` and `j`
                permuted_matrix[:, i] = m[:, j]
                permuted_matrix[:, j] = m[:, i]

                # sort by rows
                permuted_matrix[:] = permuted_matrix[argsort_2D_array(permuted_matrix)]

                # transpose since we need column-wise ordering in the next step
                permutations.append(permuted_matrix.T)

        # sort the matrices from smallest to largest
        sorted_matrix = argsort_ND_array(permutations)

        # add all matrices that have the largest possible value for `i' to list of options to check
        options = []
        for k in range(len(sorted_matrix)-1,-1,-1):
            if np.array_equal(permutations[sorted_matrix[k]][i], permutations[sorted_matrix[-1]][i]):
                options.append(permutations[sorted_matrix[k]].T)
            else:
                break

    # resulting options are all equivalent, take first
    matrix[:] = options[0]

def light_Pak_sort(matrix):
    '''
    Inplace modify the `matrix` to some ordering,
    when permutations of rows and columns (excluding
    the first) are allowed. The implementation of
    this function is described in chapter 2 of
    [Pak11]_. This function implements a lightweight
    version: In step (v), we only consider one, not
    all table copies with the maximized second column.

    .. note::
        This function may result in different
        orderings depending on the initial ordering.

    .. seealso::
        :func:`.iterative_sort`, :func:`.Pak_sort`

    :param matrix:
        2D array-like;
        The matrix to be canonicalized.

    '''
    for i in range(1,matrix.shape[1]):
        # sort all permutations of columns `i` and `j` (where `i`<=`j`) --> pick largest
        permutations = []
        for j in range(i,matrix.shape[1]):
            permuted_matrix = matrix.copy()

            # permute integration variables `i` and `j`
            permuted_matrix[:,i] = matrix[:,j]
            permuted_matrix[:,j] = matrix[:,i]

            # sort by rows
            permuted_matrix[:] = permuted_matrix[argsort_2D_array(permuted_matrix)]

            # transpose since we need column-wise ordering in the next step
            permutations.append(permuted_matrix.T)

        # find the largest `i`th column in `permutations` and keep that
        index_of_largest = argsort_ND_array(permutations)[-1]
        matrix[:] = permutations[index_of_largest].T # transpose back
