"""

Routines to split the integration between :math:`0`
and :math:`1`. This maps singularities from :math:`1`
to :math:`0`.

"""

from .common import Sector
from ..algebra import Polynomial, ExponentiatedPolynomial, refactorize
from ..misc import powerset, sympify_expression
import numpy as np
import math

# ********************** helper functions **********************

def remap_one_to_zero(polynomial, *indices):
    r'''
    Apply the transformation :math:`x \rightarrow 1 - x`
    to `polynomial` for the parameters of the given
    `indices`.

    :param polynomial:
        :class:`.Polynomial`;
        The polynomial to apply the transformation to.

    :param indices:
        arbitrarily many int;
        The indices of the ``polynomial.polysymbols``
        to apply the transformation to.

    Example:

    >>> from pySecDec.algebra import Polynomial
    >>> from pySecDec.decomposition.splitting import remap_one_to_zero
    >>> polysymbols = ['x0']
    >>> polynomial = Polynomial.from_expression('x0', polysymbols)
    >>> remap_one_to_zero(polynomial, 0)
     + (1) + (-1)*x0

    '''
    expolist, coefficients = [], []
    # go through each term of the original polynomial
    for expos, coeff in zip(polynomial.expolist, polynomial.coeffs):
        # the factors of the term that are not remapped along with the coefficient
        termpoly = Polynomial([[0 if m in indices else expos[m] for m in range(len(polynomial.polysymbols))]], [coeff], polynomial.polysymbols)
        # generate a polynomial for each remapped monomial of this term,
        # calculate the coefficients directly using binomial coefficients
        # and multiply the polynomials to the rest of the term
        for index in indices:
            n = expos[index]
            binomexpolist = np.zeros((n+1, len(polynomial.polysymbols)), dtype=np.int)
            binomexpolist[:,index] = range(n+1)
            binomcoeffs = [pow(-1,k) * math.factorial(n) // math.factorial(n-k) // math.factorial(k) for k in range(n+1)]
            binompoly = Polynomial(binomexpolist, binomcoeffs, polynomial.polysymbols)
            termpoly *= binompoly
        # keep track of all the resulting terms
        expolist.extend(termpoly.expolist)
        coefficients.extend(termpoly.coeffs)

    # generate the full remapped polynomial
    remapped_polynomial = Polynomial(expolist, coefficients, polynomial.polysymbols).simplify()

    if hasattr(polynomial, 'exponent'):
        # convert to `ExponentiatedPolynomial`
        remapped_polynomial = ExponentiatedPolynomial(
                                                         remapped_polynomial.expolist,
                                                         remapped_polynomial.coeffs,
                                                         polynomial.exponent,
                                                         polynomial.polysymbols,
                                                         copy=True
                                                     )

    return remapped_polynomial

def find_singular_sets_at_one(polynomial):
    '''
    Find all possible sets of parameters such that the `polynomial`'s
    constant term vanishes if these parameters are set to one.

    Example:

    >>> from pySecDec.algebra import Polynomial
    >>> from pySecDec.decomposition.splitting import find_singular_sets_at_one
    >>> polysymbols = ['x0', 'x1']
    >>> poly = Polynomial.from_expression('1 - 10*x0 - x1', polysymbols)
    >>> find_singular_sets_at_one(poly)
    [(1,)]

    :param polynomial:
        :class:`.Polynomial`;
        The polynomial to search in.

    '''
    # ignore parameters that do not appear
    indices_to_consider = []
    for i in range(polynomial.number_of_variables):
        if (polynomial.expolist[:,i] != 0).any():
            indices_to_consider.append(i)

    singular_sets = []
    for singular_set in powerset(indices_to_consider):
        poly_copy = polynomial.copy()
        poly_copy.expolist[:,singular_set] = 0
        poly_copy.simplify()
        if ( not poly_copy.has_constant_term() ) or ( len(poly_copy.coeffs) == 1 and poly_copy.coeffs[0] == 0 ):
            singular_sets.append(singular_set)
    return singular_sets

# **************************** split ****************************

def split(sector, seed, *indices):
    '''
    Split the integration interval :math:`[0,1]`
    for the parameters given by `indices`. The
    splitting point is fixed using `numpy's`
    random number generator.

    Return an iterator of :class:`.Sector` - the
    arising subsectors.

    :param sector:
        :class:`.Sector`;
        The sector to be split.

    :param seed;
        integer;
        The seed for the random number generator
        that is used to fix the splitting point.

    :param indices:
        arbitrarily many integers;
        The indices of the variables to be split.

    '''
    # get the splitting point
    #  --> generate ``len(indices)`` random numbers between ``1`` and ``20``,
    #      then split at ``<random>/20``
    rng = np.random.RandomState( int(seed) )
    splitting_point = [   int( rng.randint(1,20) ) / sympify_expression(20)  for idx in indices   ]

    def split_recursively(sector, indices, splitting_point):
        if not indices:
            yield sector.copy()
            return

        # We call this function recusively and pop the first index/splitting_value in each iteration
        index = indices[0]
        remaining_indices = indices[1:]
        splitting_value = splitting_point[0]
        splitting_point = splitting_point[1:]

        # split the parameter with index `index`
        #  - step1: make a copy with mapping "x --> 1 - x"
        remapped_cast = []
        remapped_other = []
        for prod in sector.cast:
            prod = prod.copy()
            mono, poly = prod.factors

            # multiply the monomial back to the polynomial and refactorize later
            poly.expolist[:,index] += mono.expolist[0,index]
            mono.expolist[:,index] = 0
            prod.factors[1] = remap_one_to_zero(poly, index)
            refactorize(prod, index)
            remapped_cast.append(prod)

        for poly in sector.other:
            remapped_other.append( remap_one_to_zero(poly, index) )

        remapped_Jacobian = remap_one_to_zero(sector.Jacobian, index)

        subsector0 = sector.copy()
        subsector1 = Sector(remapped_cast, remapped_other, remapped_Jacobian)

        #  - step2: rescale the integration variable
        def multiply_by(polynomial, number, index):
            replaced_polynomial = polynomial.replace(index, number)
            polynomial.coeffs = replaced_polynomial.coeffs

        for sector,remapping_factor in zip(
                                              [     subsector0    ,     subsector1      ],
                                              [  splitting_value  ,  1-splitting_value  ]
                                          ):
            sector.Jacobian *= remapping_factor
            multiply_by(sector.Jacobian, remapping_factor, index)
            for prod in sector.cast:
                mono, poly = prod.factors
                multiply_by(mono, remapping_factor, index)
                multiply_by(poly, remapping_factor, index)
            for poly in sector.other:
                multiply_by(poly, remapping_factor, index)

        # recursively call `split` with the `remaining_indices`
        for sector in (subsector0, subsector1):
            for subsubsector in split_recursively(sector, remaining_indices, splitting_point):
                yield subsubsector

    return split_recursively(sector, indices, splitting_point)

def split_singular(sector, seed, indices=[]):
    '''
    Split the integration interval :math:`[0,1]`
    for the parameters that can lead to singularities
    at one for the polynomials in ``sector.cast``.

    Return an iterator of :class:`.Sector` - the
    arising subsectors.

    :param sector:
        :class:`.Sector`;
        The sector to be split.

    :param seed:
        integer;
        The seed for the random number generator
        that is used to fix the splitting point.

    :param indices:
        iterables of integers;
        The indices of the variables to be split
        if required. An empty iterator means that
        all variables may potentially be split.

    '''
    # find the parameters to be split
    singular_parameters = set()
    for product in sector.cast:
        polynomial = product.factors[1]
        for singular_set in find_singular_sets_at_one(polynomial):
            singular_parameters.update(singular_set)

    # restrict splitting to selected `indices`
    if indices:
        singular_parameters.intersection_update(indices)

    # order the `singular_parameters`
    singular_parameters = list(singular_parameters)
    singular_parameters.sort()

    # split the `sector` in the `singular_parameters`
    return split(sector, seed, *singular_parameters)
