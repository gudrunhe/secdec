"""
Expansion
---------

Routines to series expand singular and nonsingular expressions.

"""

from .algebra import Product, Polynomial, ExponentiatedPolynomial
from .misc import sympify_symbols, flatten, sympify_expression
from numpy import iterable
import numpy as np
import sympy as sp

_sympy_zero = sympify_expression(0)

class OrderError(ValueError):
    '''
    This exception is raised if an expansion
    to a lower than the lowest order of an
    expression is requested.

    '''
    pass

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
    expression_variable_set_to_zero = expression.replace(index, 0)
    coeffs = [expression_variable_set_to_zero]
    inverse_i_factorial = sympify_expression(1)
    for order_i in range(order):
        inverse_i_factorial /= order_i + 1
        expression = expression.simplify().derive(index)
        coeffs.append( expression.replace(index, 0) * inverse_i_factorial )

    return Polynomial(expolist, np.array(coeffs), expression.symbols, copy=False)

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
    num_min_exp = min(numerator.expolist[:,index])
    den_min_exp = min(denominator.expolist[:,index])
    numerator.expolist[:,index] -= num_min_exp
    denominator.expolist[:,index] -= den_min_exp
    highest_pole = - num_min_exp + den_min_exp

    if order < -highest_pole:
        raise OrderError('The lowest order (%i) is higher than the requested order (%i)' %(-highest_pole,order))

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
    nonsingular_series_zero_coeffs = []

    # nonsingular expansion order is shifted by the (potentially) singular prefactor (`highest_pole`)
    i_factorial = 1
    for i in range(order + highest_pole):
        i_factorial *= i + 1

        numerator, denominator = derive_polyrational(numerator, denominator, index)
        this_order_numerator = numerator.replace(index, 0).simplify()

        # if the numerator is zero, skip calculating denominator
        # remember that this term was zero (may need to add it if we encounter a non-zero term later in the series)
        if (this_order_numerator.coeffs == 0).all():
            nonsingular_series_zero_coeffs.append(Product(this_order_numerator))
            continue

        # we encountered a non-zero term
        # if there are zero terms prior to this non-zero term, add them to the list of coeffs now
        if nonsingular_series_zero_coeffs:
            nonsingular_series_coeffs = nonsingular_series_coeffs + nonsingular_series_zero_coeffs
            nonsingular_series_zero_coeffs = []

        # process the denominator
        this_order_denominator = denominator.replace(index, 0).simplify()

        # convert denominator to ``<polynomial>**-1``
        this_order_denominator = ExponentiatedPolynomial(this_order_denominator.expolist, this_order_denominator.coeffs  * i_factorial, -1, this_order_denominator.polysymbols, copy=False)

        nonsingular_series_coeffs.append(Product(this_order_numerator, this_order_denominator))

    # if last terms are zero, adjust length of expolist
    if nonsingular_series_zero_coeffs:
        nonsingular_series_expolist  = nonsingular_series_expolist[:len(nonsingular_series_coeffs)]

    return Polynomial(nonsingular_series_expolist, nonsingular_series_coeffs, numerator.polysymbols)

def _expand_and_flatten(expression, indices, orders, expansion_one_variable):
    '''
    Expand `expression` in each variable passed via
    `indices` using the function `expansion_one_variable`.
    Flatten the resulting polynomial using the function
    :func:`pySecDec.misc.flatten`.

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
    if len(indices) == 1:
        return expansion
    else:
        return flatten(expansion, len(indices) - 1)

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

    .. seealso::
        To expand more general expressions use
        :func:`.expand_sympy`.

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

def expand_sympy(expression, variables, orders):
    '''
    Expand a sympy expression in the `variables`
    to given `orders`. Return the expansion as
    nested :class:`pySecDec.algebra.Polynomial`.

    .. seealso::
        This function is a generalization of
        :func:`.expand_singular`.

    :param expression:
        string or sympy expression;
        The expression to be expanded

    :param variables:
        iterable of strings or sympy symbols;
        The variables to expand the `expression`
        in.

    :param orders:
        iterable of integers;
        The orders to expand to.

    '''
    # basic consistency checks
    if not isinstance(expression, sp.Expr):
        expression = sympify_expression(expression)
    variables = sympify_symbols(variables, 'All `variables` must be symbols.')
    orders = np.asarray(orders)
    assert len(orders.shape) == 1, '`orders` must be vector-like'
    assert len(variables) == len(orders), 'The number of variables (%i) must equal the number of orders (%i).' % (len(variables), len(orders))

    def recursion(expression, variables, orders, index):
        variable = variables[index]
        order = orders[index]

        # get the lowest term
        expansion_generator = sp.series(expression, variable, n=None)
        lowest_order_term = next(expansion_generator)

        # find out how many more terms are required:
        #   -> step1: get the order of the lowest term
        #      try conversion to sympy polynomial --> fails for pole
        lowest_order_variable_only = (lowest_order_term/lowest_order_term.subs(variable,1)).cancel()
        try:
            lowest_order = sp.poly(lowest_order_variable_only, variable).monoms()[0][0]
        except sp.PolynomialError:
            # pole --> convert the inverse to polynomial
            lowest_order = - sp.poly(1/lowest_order_variable_only, variable).monoms()[0][0]

        #   -> step2: compute how many orders lie between the lowest and the requested order
        number_of_remaining_orders = order - lowest_order

        # if `number_of_remaining_orders` is negative, the requested order is lower than the lowest order
        if number_of_remaining_orders < 0:
            raise OrderError( 'The lowest order in `%s` (%i) is higher than the requested order (%i)' % (variable,lowest_order,order) )

        # compute the remaining terms
        coeffs = []
        current_term = lowest_order_term
        current_order = lowest_order
        truncated = True
        while current_order <= order:
            coeffs.append( (current_term/variable**(current_order)).expand() )
            current_order += 1
            try:
                current_term = next(expansion_generator)
            except StopIteration:
                # all higher orders are exactly zero
                truncated = False
                break
            current_term_variable_only = (current_term/current_term.subs(variable,1)).cancel()
            try:
                current_term_order = sp.poly(current_term_variable_only, variable).monoms()[0][0]
            except sp.PolynomialError:
                current_term_order = -sp.poly(1/current_term_variable_only, variable).monoms()[0][0]
            # if some intermediate orders are zero, append zeros to the expansion
            while current_order < current_term_order and current_order <= order:
                coeffs.append(_sympy_zero)
                current_order += 1

        coeffs = np.array(coeffs)

        # create the `Polynomial`
        expolist = np.zeros([len(coeffs),len(variables)], dtype=int)
        expolist[:,index] = np.arange(lowest_order, current_order)
        expansion = Polynomial(expolist, coeffs, variables, copy=False)
        expansion.truncated = truncated

        # recurse down to the next variable if any
        if index + 1 < len(variables):
            for i,coeff in enumerate(coeffs):
                coeffs[i] = recursion(coeff, variables, orders, index + 1)

        return expansion

    # initialize the recursion
    return recursion(expression, variables, orders, 0)
