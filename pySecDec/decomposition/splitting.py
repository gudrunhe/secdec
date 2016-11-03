"""

Routines to split the integration between :math:`0`
and :math:`1` at :math:`1/2`. This maps singularities
from :math:`1` to :math:`0`.

"""

from .common import Sector, refactorize
from ..algebra import Polynomial, ExponentiatedPolynomial
from ..misc import powerset
import numpy as np
import sympy as sp

_sympy_one_half = sp.sympify('1/2')

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
    coeffs = polynomial.coeffs
    monomials = [Polynomial.from_expression(s, polynomial.polysymbols) for s in polynomial.polysymbols]
    dummy_poly_coeffs = np.array([sp.sympify('coeffs_%i__'%i) for i in range(len(coeffs))])
    dummy_poly_polysymbols = ('monomials_%i__'%i for i in range(len(polynomial.polysymbols)))
    dummy_poly = Polynomial(polynomial.expolist.copy(), dummy_poly_coeffs,
                                        dummy_poly_polysymbols, copy=False)
    code_new_poly = str(dummy_poly).replace('__',']').replace('_', '[')
    for index in indices:
        code_new_poly = code_new_poly.replace('monomials[%i]'%index, '(1-monomials[%i])'%index)

    remapped_polynomial = eval(code_new_poly)

    # make sure that `remapped_polynomial` is of type Polynomial
    if not isinstance(remapped_polynomial, Polynomial):
        remapped_polynomial = Polynomial(np.zeros([1,len(monomials)],dtype=int), np.array([remapped_polynomial]), polynomial.polysymbols, copy=1)

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
    singular_sets = []
    for singular_set in powerset(range(len(polynomial.polysymbols))):
        poly_copy = polynomial.copy()
        poly_copy.expolist[:,singular_set] = 0
        poly_copy.simplify()
        if ( not poly_copy.has_constant_term() ) or ( len(poly_copy.coeffs) == 1 and poly_copy.coeffs[0] == 0 ):
            singular_sets.append(singular_set)
    return singular_sets

# ************************ split at 1/2 ************************

def split(sector, *indices):
    '''
    Split the integration interval :math:`[0,1]`
    at :math:`1/2` for the parameters given by
    `indices`.

    Return an iterator of :class:`.Sector` - the
    arising subsectors.

    :param sector:
        :class:`.Sector`;
        The sector to be split.

    :param indices:
        integers;
        The indices of the variables to be split.

    '''
    # Jacobian must not depend on the variable denoted by index.
    for index in indices:
        assert (sector.Jacobian.expolist[:,index] == 0).all(), "``sector.Jacobian`` must not depend on the variables denoted by `indices`."

    def split_recursively(sector, indices):
        if not indices:
            yield sector.copy()
            return

        # We call this function recusively and pop the first index in each iteration
        index = indices[0]
        remaining_indices = indices[1:]

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

        subsector0 = sector.copy()
        subsector1 = Sector(remapped_cast, remapped_other, sector.Jacobian)

        #  - step2: remap "x --> x/2"
        def divide_by_two(polynomial, index):
            replaced_polynomial = polynomial.replace(index, _sympy_one_half)
            polynomial.coeffs = replaced_polynomial.coeffs

        for sector in (subsector0, subsector1):
            sector.Jacobian *= _sympy_one_half
            for prod in sector.cast:
                mono, poly = prod.factors
                divide_by_two(mono, index)
                divide_by_two(poly, index)
            for poly in sector.other:
                divide_by_two(poly, index)

        # recursively call `split` with the `remaining_indices`
        for sector in (subsector0, subsector1):
            for subsubsector in split_recursively(sector, remaining_indices):
                yield subsubsector

    return split_recursively(sector, indices)

def split_singular(sector):
    '''
    Split the integration interval :math:`[0,1]`
    at :math:`1/2` for the parameters that can
    lead to singularities at one for the
    polynomials in ``sector.cast``.

    Return an iterator of :class:`.Sector` - the
    arising subsectors.

    :param sector:
        :class:`.Sector`;
        The sector to be split.

    '''
    # find the parameters to be split
    singular_parameters = set()
    for product in sector.cast:
        polynomial = product.factors[1]
        for singular_set in find_singular_sets_at_one(polynomial):
            singular_parameters.update(singular_set)

    # order the `singular_parameters`
    singular_parameters = list(singular_parameters)
    singular_parameters.sort()

    # split the `sector` in the `singular_parameters`
    return split(sector, *singular_parameters)
