"""miscellaneous routines"""

from itertools import chain,combinations
import sympy as sp
import numpy as np

def powerset(iterable, exclude_empty=False, stride=1):
    """
    Return an iterator over the powerset of a given set.
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    :param iterable:
        iterable;
        The set to generate the powerset for.

    :param exclude_empty:
        bool, optional;
        If True, skip the empty set in the powerset.
        Default is False.

    :param stride:
        integer;
        Only generate sets that have a multiple of
        `stride` elements.
        powerset([1,2,3], stride=2) --> () (1,2) (1,3) (2,3)

    """
    # taken from python's own documentation
    s = list(iterable)
    powerset_iterator = iter(chain.from_iterable(combinations(s, r) for r in range(0,len(s)+1,stride)))
    if exclude_empty:
        # The first element of the iterator is the empty set -> discard
        next(powerset_iterator)
    return powerset_iterator

def missing(full, part):
    '''
    Return the elements of `full` that have to
    be appended to `part` such that `full` is
    a subset of `part`.
    missing([1,2,3], [1]) --> [2,3]
    missing([1,2,3], [1,'a']) --> [2,3]

    :param full:
        iterable;
        The set of elements to complete `part`
        with.

    :param part:
        iterable;
        The set to be completed to a superset
        of `full`.

    '''
    part = list(part)
    missing = []
    for item in full:
        if item not in part:
            missing.append(item)
    return missing

def det(M):
    '''
    Calculate the determinant of a matrix.

    :param M:
        a square-matrix-like array;

    '''
    M = np.asarray(M)
    assert len(M.shape) == 2, "`M` must be two dimensional"
    assert M.shape[0] == M.shape[1], "`M` must be a square matrix"
    D = M.shape[0]

    # Use sympy to calculate the determinant of a generic DxD Matrix, e.g.
    # [[M_0_0__, M_0_1__]
    #  [M_1_0__, M_1_1__]]
    generic_m = sp.Matrix([["M_%i_%i__" % (i,j) for j in range(D)] for i in range(D)])
    generic_det = generic_m.det(method='berkowitz').expand()

    # convert sympy output to python executable code; i.e. M_i_j__ --> M[i,j]
    algebraic_det = str(generic_det).replace('M_','M[').replace('__',']').replace('_',',')

    # execute the expression and return the result
    return eval(algebraic_det)

def adjugate(M):
    '''
    Calculate the adjugate of a matrix.

    :param M:
         a square-matrix-like array;

    '''
    M = np.asarray(M)
    assert len(M.shape) == 2, "`M` must be two dimensional"
    assert M.shape[0] == M.shape[1], "`M` must be a square matrix"
    D = M.shape[0]

    if D == 1:
        # whatever the entry of a 1x1 matrix is, its adjugate is [[1]]
        return np.array([[1]], dtype=M.dtype)

    # Use sympy to calculate the adjugate of a generic DxD Matrix, e.g.
    # [[M_0_0__, M_0_1__]
    #  [M_1_0__, M_1_1__]]
    generic_m = sp.Matrix([["M_%i_%i__" %(i,j) for j in range(D)] for i in range(D)])
    generic_adjugate = generic_m.adjugate().expand()

    # convert sympy output to python executable code; i.e. M_i_j__ --> M[i,j]
    algebraic_adjugate = [[str(generic_adjugate[i,j]).replace('M_','M[').replace('__',']').replace('_',',') for j in range(D)] for i in range(D)]

    # generate adjugate of M
    adjugate_M = np.empty((D,D), dtype=M.dtype)
    for i in range(D):
        for j in range(D):
            adjugate_M[i,j] = eval(algebraic_adjugate[i][j])

    return adjugate_M

def argsort_2D_array(array):
    '''
    Sort a 2D array according to its row entries.
    The idea is to bring identical rows together.

    Example:
        [[1,2,3],         [[1,2,3],
         [2,3,4],  --->    [1,2,3],
         [1,2,3]]          [2,3,4]]

    Return the indices like numpy's :func:`argsort`
    would.

    :param array:
        2D array;
        The array to be argsorted.

    '''
    # create a new view
    # the "[:]" is essential; if it was missing, we could overwrite the input
    array = np.asarray(array)[:]

    assert len(array.shape) == 2, "`array` must be two dimensional"

    # reinterpret each column as a single data type (C struct)
    # see also `record arrays` in the numpy documentation
    array.dtype = [('column' + str(i), array.dtype) for i in range(array.shape[1])]

    array = array.flatten()

    return np.argsort(array, kind='mergesort')
