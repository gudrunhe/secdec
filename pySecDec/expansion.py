"Routines to series expand singular and nonsingular expressions"

from .algebra import Product, Sum, Polynomial, ExponentiatedPolynomial, replace, _get_symbols
from numpy import iterable
import numpy as np
import sympy as sp

# ------------------------ private functions ------------------------

def _expand_Taylor_step(expression, index, order):
    r'''
    Series/Taylor expand a nonsingular `expression` around
    zero.

    :param expression:
        an expression composed of the types defined in
        the module :mod:`.algebra`;
        The expression to be series expanded.

    :param index:
        integer;
        The index of the parameter to expand.

    :param order:
        nonnegative integer;
        The order to which the expansion is to be calculated.

    '''
    int_order = int(order)
    assert int_order == order and int_order >= 0, "`order` must be a nonnegative integer"
    order = int_order
    N = expression.number_of_variables

    # Construct expolist of the Tatlor polynomial
    expolist = np.zeros((1 + order, N), dtype=int)
    expolist[:,index] = np.arange(1 + order)

    # Construct coefficients of the Taylor polynomial
    expression_variable_set_to_zero = replace(expression, index, 0).simplify()
    coeffs = [expression_variable_set_to_zero]
    for order_i in range(order):
        expression = expression.derive(index).simplify()
        coeffs.append( replace(expression, index, 0).simplify() )

    return Polynomial(expolist, coeffs, _get_symbols(expression))

def _expand_singular_step(product, index, order):
    r'''
    Series expand a potentially singular expression of the form

    .. math::
        {
            \frac{a_N \epsilon_0 + b_N \epsilon_1 + ...}
                 {a_D \epsilon_0 + b_D \epsilon_1 + ...}
        }

    Return a :class:`.algebra.Polynomial` with coefficients
    of :class:`.algebra.Product` of the same form
    as above - the series expansion.

    :param product:
        :class:`.algebra.Product` of the form
        ``<numerator polynomial> * <denominator polynomial> ** -1``;
        The expression to be series expanded.

    :param index:
        integer;
        The index of the parameter to expand in.

    :param order:
        integer;
        The order to which the expansion is to be calculated.

    '''
    # must have a rational polynomial (polynomial product of the form p * p**-1) in the first arg
    if type(product) is not Product:
        raise TypeError('`product` must be a `Product`')
    if len(product.factors) != 2:
        raise TypeError('`product` must consist of exactly two factors')
    if type(product.factors[0]) is not Polynomial:
        raise TypeError('The first factor of `product` must be a `Polynomial` (not a subtype)')
    if not isinstance(product.factors[1], ExponentiatedPolynomial):
        raise TypeError('The second factor of `product` must be an `ExponentiatedPolynomial` with ``exponent=-1``')
    if product.factors[1].exponent != -1:
        raise TypeError('The second factor of `product` must be an `ExponentiatedPolynomial` with ``exponent=-1``')

    N = product.number_of_variables
    numerator = product.factors[0].copy()
    denominator = Polynomial(product.factors[1].expolist, product.factors[1].coeffs, product.factors[1].polysymbols) # convert to `Polynomial` (without exponent)

    # factorize overall epsilon from numerator and denominator such that it is nonzero for epsilon -> 0
    #    => can do ordinary Taylor expansion on the non singular factor
    highest_pole = 0
    while denominator.becomes_zero_for([index]):
        highest_pole += 1
        denominator.expolist[:,index] -= 1
    while numerator.becomes_zero_for([index]):
        highest_pole -= 1
        numerator.expolist[:,index] -= 1

    if order < -highest_pole:
        raise ValueError('The lowest order (%i) is higher than the requested order (%i)' %(-highest_pole,order))

    def derive_polyrational(numerator, denominator, index):
        '''
        Return the derivative of a polyrational function
        as polyrational function (numerator, denominator).

        '''
        new_numerator = numerator.derive(index) * denominator - denominator.derive(index) * numerator
        new_denominator = denominator * denominator
        return new_numerator, new_denominator


    # Taylor expansion of the nonsingular rational polynomial defined by ``numerator / denominator``
    this_order_numerator = replace(numerator, index, 0)
    this_order_denominator = replace(denominator, index, 0)

    # convert denominator to ``<polynomial>**-1``
    this_order_denominator = ExponentiatedPolynomial(this_order_denominator.expolist, this_order_denominator.coeffs, -1, this_order_denominator.polysymbols)

    nonsingular_series_coeffs = [Product(this_order_numerator, this_order_denominator)]
    nonsingular_series_expolist = np.zeros((1 + order + highest_pole, N), dtype=int)
    nonsingular_series_expolist[:,index] = np.arange(1 + order + highest_pole) - highest_pole

    # nonsingular expansion order is shifted by the (potentially) singular prefactor (`highest_pole`)
    for i in range(order + highest_pole):
        numerator, denominator = derive_polyrational(numerator, denominator, index)
        this_order_numerator = replace(numerator, index, 0)

        # if the numerator is zero at one order, all higher orders vanish as well
        if (this_order_numerator.coeffs == 0).all():
            nonsingular_series_expolist = nonsingular_series_expolist[:i+1]
            break

        this_order_denominator = replace(denominator, index, 0)

        # convert denominator to ``<polynomial>**-1``
        this_order_denominator = ExponentiatedPolynomial(this_order_denominator.expolist, this_order_denominator.coeffs, -1, this_order_denominator.polysymbols)

        nonsingular_series_coeffs.append(Product(this_order_numerator, this_order_denominator))

    return Polynomial(nonsingular_series_expolist, nonsingular_series_coeffs, numerator.polysymbols)

def _flatten(polynomial):
    '''
    Convert the output of :func:`_expand_singular_step`
    to a polynomial whose coefficients only consist of
    numbers or symbols.

    :param polynomial:
        :class:`pySecDec.algebra.Polynomial`;
        The polynomial to "flatten".

    '''
    def recursively_sympify(expr):
        assert isinstance(expr, Polynomial)
        for i,coeff in enumerate(expr.coeffs):
            if not isinstance(expr.coeffs[i], Polynomial):
                expr.coeffs[i] = sp.sympify(coeff)
            else:
                recursively_sympify(coeff)

    def flatten_recursively(expr):
        '''
        Expand `Polynomial`s in the coefficients into
        the top polynomial

        '''
        assert isinstance(expr, Polynomial)
        outpoly = 0
        for i,(exponents,coeff) in enumerate(zip(expr.expolist,expr.coeffs)):
            if isinstance(coeff, Polynomial):
                coeff = flatten_recursively(coeff)
            monomial = Polynomial([exponents],[1], expr.polysymbols)
            outpoly += monomial * coeff
        return outpoly

    recursively_sympify(polynomial)
    return flatten_recursively(polynomial)

def _expand_and_flatten(expression, indices, orders, expansion_one_variable):
    '''
    Expand `expression` in each variable passed via
    `indices` using the function `expansion_one_variable`.
    Flatten the resulting polynomial using the function
    :func:`_flatten`.

    '''
    # basic consistency check
    # Note: More checks are done by `_expand_Taylor_step` or `_expand_singular_step`
    indices = list(indices) if iterable(indices) else [indices]
    orders = list(orders) if iterable(orders) else [orders]
    assert len(indices) == len(orders), '`indices` and `orders` must have the same length'

    # reverse lists because popping the last element of a list is cheaper (order one) than the first (order N)
    indices.reverse()
    orders.reverse()

    def recursive_expansion(expression, indices, orders):
        # make a copy
        indices, orders = list(indices), list(orders)

        # pop the last index and order (lists have been reversed) and expand in that parameter
        index, order = indices.pop(), orders.pop()
        expansion = expansion_one_variable(expression, index, order)

        if indices:
            assert orders
            # recursively expand the other parameters in the coefficients
            for i, coeff in enumerate(expansion.coeffs):
                expansion.coeffs[i] = recursive_expansion(coeff, indices, orders)

        return expansion


    expansion = recursive_expansion(expression, indices, orders)
    return _flatten(expansion)

# -------------------- end of private functions --------------------

def expand_singular(product, indices, orders):
    r'''
    Series expand a potentially singular expression of the form

    .. math::
        {
            \frac{a_N \epsilon_0 + b_N \epsilon_1 + ...}
                 {a_D \epsilon_0 + b_D \epsilon_1 + ...}
        }

    Return a :class:`.algebra.Polynomial` - the series expansion.

    :param product:
        :class:`.algebra.Product` of the form
        ``<numerator polynomial> * <denominator polynomial> ** -1``;
        The expression to be series expanded.

    :param indices:
        integer or iterable of integers;
        The indices of the parameters to expand. The ordering of
        the indices defines the ordering of the expansion.

    :param order:
        integer or iterable of integers;
        The order to which the expansion is to be calculated.

    '''
    return _expand_and_flatten(product, indices, orders, _expand_singular_step)

def expand_Taylor(expression, indices, orders):
    r'''
    Series/Taylor expand a nonsingular `expression` around
    zero.

    Return a :class:`.algebra.Polynomial` - the series expansion.

    :param expression:
        an expression composed of the types defined in
        the module :mod:`.algebra`;
        The expression to be series expanded.

    :param indices:
        integer or iterable of integers;
        The indices of the parameters to expand. The ordering of
        the indices defines the ordering of the expansion.

    :param order:
        integer or iterable of integers;
        The order to which the expansion is to be calculated.

    '''
    return _expand_and_flatten(expression, indices, orders, _expand_Taylor_step)
