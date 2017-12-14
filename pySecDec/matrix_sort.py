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

def Pak_sort(matrix, *indices):
    '''
    Inplace modify the `matrix` to some canonical ordering,
    when permutations of rows and columns are allowed.

    The `indices` parameter can contain a list of lists
    of column indices. Only the columns present in the
    same list are swapped with each other.

    The implementation of this function is described
    in chapter 2 of [Pak11]_.

    .. note::
        If not all indices are considered the resulting
        matrix may not be canonical.

    .. seealso::
        :func:`.iterative_sort`, :func:`.light_Pak_sort`

    :param matrix:
        2D array-like;
        The matrix to be canonicalized.

    :param indices:
        arbitrarily many iterables of non-negative integers;
        The groups of columns to permute.
        Default: ``range(1,matrix.shape[1])``

    '''
    if indices:
        # Always consider the smallest indices first
        indices = map(sorted,indices)
    else:
        indices = [range(1,matrix.shape[1])]

    options = [matrix]
    for set_of_indices in indices:
        for idx, i in enumerate(set_of_indices):
            remaining_indices = set_of_indices[idx:]
            permutations = []
            for m in options:
                for j in remaining_indices:
                    permuted_matrix = m.copy()

                    # permute integration variables `column_to_swap` and `j`
                    permuted_matrix[:, i] = m[:, j]
                    permuted_matrix[:, j] = m[:, i]

                    # sort by rows
                    permuted_matrix = permuted_matrix[argsort_2D_array(permuted_matrix)]

                    # transpose since we need column-wise ordering in the next step
                    permutations.append(permuted_matrix.T)

            # sort the matrices from smallest to largest
            sorted_matrix = argsort_ND_array(permutations)

            # add all matrices that have the smallest possible value for `i' to list of options to check
            options = []
            for k in range(len(sorted_matrix)):
                if np.array_equal(permutations[sorted_matrix[k]][i], permutations[sorted_matrix[0]][i]):
                    for mat in options:
                        if np.array_equal(mat, permutations[sorted_matrix[k]].T):
                            break
                    else:
                        options.append(permutations[sorted_matrix[k]].T)
                else:
                    break

    # Pick one possible option
    matrix[:] = options[0]

def light_Pak_sort(matrix):
    '''
    Inplace modify the `matrix` to some ordering,
    when permutations of rows and columns (excluding
    the first) are allowed. The implementation of
    this function is described in chapter 2 of
    [Pak11]_. This function implements a lightweight
    version: In step (v), we only consider one, not
    all table copies with the minimized second column.

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

        # find the smallest `i`th column in `permutations` and keep that
        index_of_smallest = argsort_ND_array(permutations)[0]
        matrix[:] = permutations[index_of_smallest].T # transpose back
