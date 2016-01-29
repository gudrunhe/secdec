"Routines to series expand singular and nonsingular expressions"

from .algebra import Product, Sum, Polynomial, ExponentiatedPolynomial
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

    # Construct expolist of the Taylor polynomial
    expolist = np.zeros((1 + order, N), dtype=int)
    expolist[:,index] = np.arange(1 + order)

    # Construct coefficients of the Taylor polynomial
    expression_variable_set_to_zero = expression.replace(index, 0).simplify()
    coeffs = [expression_variable_set_to_zero]
    for order_i in range(order):
        expression = expression.derive(index).simplify()
        coeffs.append( expression.replace(index, 0).simplify() )

    return Polynomial(expolist, coeffs, expression.symbols)

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
        :class:`.algebra.Product` with factors of the form
        ``<polynomial>`` or ``<polynomial> ** -1``;
        The expression to be series expanded.

    :param index:
        integer;
        The index of the parameter to expand in.

    :param order:
        integer;
        The order to which the expansion is to be calculated.

    '''
    N = product.number_of_variables
    symbols = product.symbols
    numerator = Polynomial(np.zeros([1,N], dtype=int), np.array([1]), symbols, copy=False)
    denominator = Polynomial(np.zeros([1,N], dtype=int), np.array([1]), symbols, copy=False)

    # must have a rational polynomial (product with factors of the form <p> and <p**-1>)
    if type(product) is not Product:
        raise TypeError('`product` must be a `Product`')
    for factor in product.factors:
        # type `Polynomial` --> numerator
        if type(factor) is Polynomial:
            numerator *= factor
        # type `ExponentiatedPolynomial` with ``exponent==-1`` --> denominator
        elif type(factor) is ExponentiatedPolynomial:
            if factor.exponent != -1:
                raise TypeError('All `factors` of `product` of type `ExponentiatedPolynomial` must fulfill ``(exponent==-1) is True``')
            denominator *= Polynomial(factor.expolist, factor.coeffs, factor.polysymbols, copy=False)
        # other type --> wtf??
        else:
            raise TypeError('All `factors` of `product` must be of type `Polynomial` or `ExponentiatedPolynomial`')

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
    this_order_numerator = numerator.replace(index, 0).simplify()
    this_order_denominator = denominator.replace(index, 0).simplify()

    # convert denominator to ``<polynomial>**-1``
    this_order_denominator = ExponentiatedPolynomial(this_order_denominator.expolist, this_order_denominator.coeffs, -1, this_order_denominator.polysymbols, copy=False)

    nonsingular_series_coeffs = [Product(this_order_numerator, this_order_denominator)]
    nonsingular_series_expolist = np.zeros((1 + order + highest_pole, N), dtype=int)
    nonsingular_series_expolist[:,index] = np.arange(1 + order + highest_pole) - highest_pole

    # nonsingular expansion order is shifted by the (potentially) singular prefactor (`highest_pole`)
    for i in range(order + highest_pole):
        numerator, denominator = derive_polyrational(numerator, denominator, index)
        this_order_numerator = numerator.replace(index, 0).simplify()

        # if the numerator is zero at one order, all higher orders vanish as well
        if (this_order_numerator.coeffs == 0).all():
            nonsingular_series_expolist = nonsingular_series_expolist[:i+1]
            break

        this_order_denominator = denominator.replace(index, 0).simplify()

        # convert denominator to ``<polynomial>**-1``
        this_order_denominator = ExponentiatedPolynomial(this_order_denominator.expolist, this_order_denominator.coeffs, -1, this_order_denominator.polysymbols)

        nonsingular_series_coeffs.append(Product(this_order_numerator, this_order_denominator))

    return Polynomial(nonsingular_series_expolist, nonsingular_series_coeffs, numerator.polysymbols)

def _flatten(polynomial, depth):
    '''
    Convert the output of :func:`_expand_singular_step`
    to a polynomial in the exapnsion variables.

    :param polynomial:
        :class:`pySecDec.algebra.Polynomial`;
        The polynomial to "flatten".

    :param depth:
        integer;
        The number of recursion steps.

    '''
    assert isinstance(polynomial, Polynomial)
    outpoly = 0
    if depth == 1:
        for i,(exponents,coeff) in enumerate(zip(polynomial.expolist,polynomial.coeffs)):
            # no further recursion
            monomial = Polynomial([exponents], [coeff], polynomial.polysymbols)
            outpoly += monomial
        return outpoly
    else:
        for i,(exponents,coeff) in enumerate(zip(polynomial.expolist,polynomial.coeffs)):
            coeff = _flatten(coeff, depth-1)
            monomial = Polynomial([exponents], [coeff], polynomial.polysymbols)
            outpoly += monomial
        return outpoly

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
    return _flatten(expansion, len(indices))

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
        :class:`.algebra.Product` with factors of the form
        ``<polynomial>`` and ``<polynomial> ** -1``;
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
