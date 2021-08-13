"""
Miscellaneous
-------------

Collection of general-purpose helper functions.

"""

from itertools import chain, combinations, product
import sympy as sp
import numpy as np
import warnings

def make_cpp_list(python_list):
    '''
    Convert a python list to a string to be used
    in c++ initializer list.

    Example: ``['a', 'b', 'c'] --> '"a","b","c"'``

    '''
    joined_inner_part = '","'.join(str(item) for item in python_list)
    if len(joined_inner_part) == 0:
        return ''
    else:
        return '"' + joined_inner_part + '"'

def powerset(iterable, min_length=0, stride=1):
    """
    Return an iterator over the powerset of a given set.
    ``powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3)
    (2,3) (1,2,3)``

    :param iterable:
        iterable;
        The set to generate the powerset for.

    :param min_length:
        integer, optional;
        Only generate sets with minimal given length.
        Default: ``0``.

    :param stride:
        integer;
        Only generate sets that have a multiple of
        `stride` elements.
        ``powerset([1,2,3], stride=2) --> () (1,2) (1,3)
        (2,3)``

    """
    # taken from python's own documentation
    s = list(iterable)
    powerset_iterator = iter(chain.from_iterable(combinations(s, r) for r in range(min_length,len(s)+1,stride)))
    return powerset_iterator

def rangecomb(low, high):
    '''
    Return an iterator over the occuring orders in a
    multivariate series expansion between `low` and
    `high`.

    :param low:
        vector-like array; The lowest orders.

    :param high:
        vector-like array; The highest orders.

    Example:

    >>> from pySecDec.misc import rangecomb
    >>> all_orders = rangecomb([-1,-2], [0,0])
    >>> list(all_orders)
    [(-1, -2), (-1, -1), (-1, 0), (0, -2), (0, -1), (0, 0)]

    '''
    low = np.asarray(low)
    high = np.asarray(high)

    assert len(low.shape) == 1, '``low`` is not vector like'
    assert len(high.shape) == 1, '``high`` is not vector like'
    assert len(low) == len(high), '``low`` (length %i) and ``high`` (length %i) must have the same length.' % (len(low),len(high))

    return product( *(np.arange(l,h) for l,h in zip(low,high+1)) )

def missing(full, part):
    '''
    Return the elements in `full` that are
    not contained in `part`. Raise `ValueError`
    if an element is in `part` but not in `full`.
    ``missing([1,2,3], [1]) --> [2,3]``
    ``missing([1,2,3,1], [1,2]) --> [3,1]``
    ``missing([1,2,3], [1,'a']) --> ValueError``

    :param full:
        iterable;
        The set of elements to complete `part`
        with.

    :param part:
        iterable;
        The set to be completed to a superset
        of `full`.

    '''
    missing = list(full)
    for item in part:
        missing.remove(item)
    return missing

def all_pairs(iterable):
    '''
    Return all possible pairs of a given set.
    ``all_pairs([1,2,3,4]) --> [(1,2),(3,4)]
    [(1,3),(2,4)] [(1,4),(2,3)]``

    :param iterable:
        iterable;
        The set to be split into all possible pairs.

    '''
    # the following iterative routine is taken from
    # http://stackoverflow.com/questions/5360220/how-to-split-a-list-into-pairs-in-all-possible-ways
    def all_pairs_recursion(lst):
        if len(lst) < 2:
            yield lst
            return
        a = lst[0]
        for i in range(1,len(lst)):
            pair = (a,lst[i])
            for rest in all_pairs_recursion(lst[1:i]+lst[i+1:]):
                yield [pair] + rest

    lst = iterable if isinstance(iterable, list) else list(iterable)
    assert len(lst) % 2 == 0, '`iterable` must have even length'
    return all_pairs_recursion(lst)

def parallel_det(M, pool):
    '''
    Calculate the determinant of a matrix in parallel.

    :param M:
        a square-matrix-like array;

    :param pool:
        :class:`multiprocessing.Pool`;
        The pool to be used.

    Example:

    >>> from pySecDec.misc import parallel_det
    >>> from multiprocessing import Pool
    >>> from sympy import sympify
    >>> M = [['m11','m12','m13','m14'],
    ...      ['m21','m22','m23','m24'],
    ...      ['m31','m32','m33','m34'],
    ...      ['m41','m42','m43','m44']]
    >>> M = sympify(M)
    >>> parallel_det(M, Pool(2)) # 2 processes
    m11*(m22*(m33*m44 - m34*m43) - m23*(m32*m44 - m34*m42) + m24*(m32*m43 - m33*m42)) - m12*(m21*(m33*m44 - m34*m43) - m23*(m31*m44 - m34*m41) + m24*(m31*m43 - m33*m41)) + m13*(m21*(m32*m44 - m34*m42) - m22*(m31*m44 - m34*m41) + m24*(m31*m42 - m32*m41)) - m14*(m21*(m32*m43 - m33*m42) - m22*(m31*m43 - m33*m41) + m23*(m31*m42 - m32*m41))

    '''
    M = np.asarray(M)
    assert len(M.shape) == 2, "`M` must be two dimensional"
    assert M.shape[0] == M.shape[1], "`M` must be a square matrix"
    D = M.shape[0]

    # stopping criterion for recursion
    if D == 1:
        return M[0,0]

    # fast check if an integer is even
    is_even = lambda x: x == (x >> 1 << 1)

    def sub_indices(i):
        'Return ``list(range(D))`` omitting `i`'
        sub_indices = list(range(D))
        sub_indices.remove(i)
        return sub_indices

    # resolve the first row of the matix
    sub_Ms = (M[[[k] for k in range(1,D)],sub_indices(i)] for i in range(D))
    sub_dets = pool.imap(det, sub_Ms)
    result = 0
    for i,subdet in enumerate(sub_dets):
        term = M[0,i] * subdet
        if is_even(i):
            result += term
        else:
            result -= term
    return result

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

    # stopping criterion for recursion
    if D == 1:
        return M[0,0]

    # fast check if an integer is even
    is_even = lambda x: x == (x >> 1 << 1)

    def sub_indices(i):
        'Return ``list(range(D))`` omitting `i`'
        sub_indices = list(range(D))
        sub_indices.remove(i)
        return sub_indices

    # resolve the first row of the matix
    for i in range(D):
        sub_M = M[[[k] for k in range(1,D)],sub_indices(i)]
        term = M[0,i] * det(sub_M) # recursion
        if i == 0:
            result = term
        elif is_even(i):
            result += term
        else:
            result -= term
    return result

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
    r'''
    Sort a 2D array according to its row entries.
    The idea is to bring identical rows together.

    .. seealso::
        If your array is not two dimesional
        use :func:`.argsort_ND_array`.

    Example:
        +-------+--------+-------+
        | input |        |sorted |
        +=======+========+=======+
        | 1 2 3 |        | 1 2 3 |
        +-------+--------+-------+
        | 2 3 4 |        | 1 2 3 |
        +-------+--------+-------+
        | 1 2 3 |        | 2 3 4 |
        +-------+--------+-------+

    Return the indices like numpy's :func:`argsort`
    would.

    :param array:
        2D array;
        The array to be argsorted.

    '''
    # create a new view
    # Copying the array is essential, otherwise we could overwrite the input.
    # Moreover, sometimes the algorithm crashed if the array is not contiguous
    array = np.asarray(array).copy()

    assert len(array.shape) == 2, "`array` must be two dimensional"

    # reinterpret each column as a single data type (C struct)
    # see also `record arrays` in the numpy documentation
    array.dtype = [('column' + str(i), array.dtype) for i in range(array.shape[1])]

    array = array.flatten()

    return np.argsort(array, kind='mergesort')

def argsort_ND_array(array):
    '''
    Like :func:`.argsort_2D_array`, this
    function groups identical entries in
    an array with any dimensionality greater
    than (or equal to) two together.

    Return the indices like numpy's :func:`argsort`
    would.

    .. seealso::
        :func:`.argsort_2D_array`

    :param array:
        ND array, :math:`N>=2`;
        The array to be argsorted.

    '''
    # create a new view
    # the "[:]" is essential; if it was missing, we could overwrite the input
    array = np.asarray(array)[:]
    assert len(array.shape) >= 2, "`array` must be two or higher dimensional"
    shape_2D = (array.shape[0], np.prod(array.shape[1:]))
    array_2D = array.reshape(shape_2D)
    return argsort_2D_array(array_2D)

def cached_property(method):
    '''
    Like the builtin `property` to be used as decorator
    but the method is only called once per instance.

    Example:

    .. code-block:: python

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

    '''
    def wrapped_method(obj):
        try:
            # try to return the result from cache
            return obj.__dict__[method.__name__]
        except KeyError:
            # call the function once and for all
            result = method(obj)
            obj.__dict__[method.__name__] = result
            return result
    # make the method a property
    return property(wrapped_method)

def doc(docstring):
    '''
    Decorator that replaces a function's docstring
    with `docstring`.

    Example:

    .. code-block:: python

        @doc('documentation of `some_funcion`')
        def some_function(*args, **kwargs):
            pass

    '''
    def add_doc(function):
        function.__doc__ = docstring
        return function
    return add_doc

def flatten(polynomial, depth=np.inf):
    '''
    Convert nested polynomials; i.e. polynomials
    that have polynomials in their coefficients
    to one single polynomial.

    :param polynomial:
        :class:`pySecDec.algebra.Polynomial`;
        The polynomial to "flatten".

    :param depth:
        integer;
        The maximum number of recursion steps.
        If not provided, stop if the coefficient
        is not a :class:`pySecDec.algebra.Polynomial`.

    '''
    from .algebra import Polynomial
    assert isinstance(polynomial, Polynomial)
    all_exponents = []
    all_coeffs = []
    for i,(exponents,coeff) in enumerate(zip(polynomial.expolist,polynomial.coeffs)):
        if type(coeff) is Polynomial: # stopping criterion for recursion
            if depth > 1: # stopping criterion for recursion
                coeff = flatten(coeff, depth-1)
            all_coeffs.extend(coeff.coeffs)
            all_exponents.append(exponents + coeff.expolist)
        else:
            all_coeffs.append(coeff)
            all_exponents.append(exponents)
    return Polynomial(np.vstack(all_exponents), np.array(all_coeffs), polynomial.polysymbols, copy=False)

def sympify_expression(a):
    '''
    A helper function for converting objects to
    sympy expressions.

    First try to convert object to a sympy expression,
    if this fails, then try to convert str(object)
    to a sympy expression

    :param a:
    The object to be converted to a sympy expression.

    :return:
    A sympy expression representing the object.
    '''
    with warnings.catch_warnings():

        # Promote Sympy "String fallback in sympify" (issue 18066) deprecation warning to error
        warnings.filterwarnings(action='error', category=sp.utilities.exceptions.SymPyDeprecationWarning, module='sympy.core.sympify')

        # Catch SymPyDeprecationWarning/SympifyError from sympify(e) and try sympify(str(e)) instead
        try:
            return sp.sympify(a)
        except (sp.utilities.exceptions.SymPyDeprecationWarning, sp.SympifyError):
            return sp.sympify(str(a))


def sympify_symbols(iterable, error_message, allow_number=False):
    '''
    `sympify` each item in `iterable` and assert
    that it is a `symbol`.

    '''
    symbols = []
    for expression in iterable:
        expression = sympify_expression(expression)
        assert expression.is_Symbol or expression.is_Number if allow_number else expression.is_Symbol, error_message
        symbols.append(expression)
    return symbols

def assert_degree_at_most_max_degree(expression, variables, max_degree, error_message):
    '''
    Assert that `expression` is a polynomial of
    degree less or equal `max_degree` in the `variables`.

    '''
    from .algebra import Polynomial
    poly = Polynomial.from_expression(expression, variables)
    assert (poly.expolist.sum(axis=1) <= max_degree).all(), error_message

def lowest_order(expression, variable):
    '''
    Find the lowest order of `expression`'s series
    expansion in `variable`.

    Example:

    >>> from pySecDec.misc import lowest_order
    >>> lowest_order('exp(eps)', 'eps')
    0
    >>> lowest_order('gamma(eps)', 'eps')
    -1

    :param expression:
        string or sympy expression;
        The expression to compute the lowest
        expansion order of.

    :param variable:
        string or sympy expression;
        The variable in which to expand.

    '''
    # convert to sympy if neccessary
    variable = sympify_symbols([variable], '`variable` must be a symbol')[0]
    if not isinstance(expression, sp.Expr):
        expression = sympify_expression(expression)

    # get the lowest term in the expansion
    lowest_expansion_term = next( sp.series(expression, variable, n=None) )

    # try conversion to sympy polynomial --> fails for pole
    try:
        lowest_expansion_term = sp.poly(lowest_expansion_term, variable)
        return lowest_expansion_term.monoms()[0][0]
    except sp.PolynomialError:
        # pole --> convert the inverse to polynomial
        highest_pole = sp.poly(lowest_expansion_term**-1, variable)
        return - highest_pole.monoms()[0][0]

def chunks(lst, n):
    '''
    Yield successive n-sized chunks from lst.

    :param lst:
        list;
        The list from which to produce chunks.

    :param n:
        integer;
        The size of the chunks to produce.

    :return:
        A list of at most length n.

    '''
    if n == 0:
        n = n + 1
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def rec_subs(expr, repl, n = 200):
    '''
    Substitute repl in expr and expand until
    expr stops changing or depth n is reached.

    :param expr:
        sympy expression;
        Expression to which the replacement
        rules are applied.

    :param repl:
        list;
        List of replacement rules.

    :param n:
        integer;
        Maximal number of repl applications.

    :return:
        Expression after substitutions.

    '''
    i = 0
    while i < n:
        expr_new = expr.expand().subs(repl)
        if expr_new-expr == 0:
            return expr_new
        else:
            expr = expr_new
            i = i + 1
    raise ValueError('Possible infinite recursion in replacement_rules (depth '+str(n)+')')
