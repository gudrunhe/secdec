"""
This module implements the main program - the function
:func:`.make_package`.

"""

from ..misc import sympify_symbols
from ..algebra import _Expression, Expression, Polynomial, ExponentiatedPolynomial
import numpy as np
import sympy as sp

# The only public object this module provides is the function `make_package`.
# The module is organized in multiple sections dedicated to specific tasks
# to be addressed while writing the c++ package.
# To find the main function `make_package`, it is easiest to full-text search
# for "def make_package(".

# ----------------------------------- parse input -----------------------------------
def _parse_expressions(expressions, polysymbols, target_type, name_of_make_argument_being_parsed):
    '''
    Convert `polynomials` from string or sympy expression
    to :class:`_Expression` or :class:`ExponentiatedPolynomial`.
    If a polynomial in `polynomials` is already of that type,
    check that the `polysymbols` are correct.
    Return the parsed polynomials and their names.

    Note: `target_type` may only be :class:`_Expression` or
          :class:`ExponentiatedPolynomial`.

    '''
    parsed = []
    for expression in expressions:
        if isinstance(expression, _Expression):
            if expression.symbols != polysymbols:
                raise ValueError('"%s" (in `%s`) depends on the wrong `symbols` (is: %s, should be: %s). Try passing it as string or sympy expression.' \
                                 % (expression, name_of_make_argument_being_parsed, expression.symbols, polysymbols))
        if target_type == ExponentiatedPolynomial:
            if type(expression) == ExponentiatedPolynomial:
                if isinstance(expression.exponent, _Expression):
                    if expression.exponent.symbols != polysymbols:
                        raise ValueError('The exponent of "%s" (in `%s`) depends on the wrong `symbols` (is: %s, should be: %s). Try passing it as string or sympy expression.' \
                                         % (expression, name_of_make_argument_being_parsed, expression.symbols, polysymbols))
                    if type(expression.exponent) != Polynomial:
                        raise ValueError('The exponent of "%s" (in `%s`) should be of type `Polynomial` but is of type `%s`. Try passing it as string or sympy expression.' \
                                          % (expression,name_of_make_argument_being_parsed,type(expression)))
                else:
                    expression.exponent = Polynomial.from_expression(expression.exponent, polysymbols)
            elif type(expression) == Polynomial:
                expression = ExponentiatedPolynomial(expression.expolist, expression.coeffs,
                                                     Polynomial.from_expression(1, polysymbols), # exponent
                                                     polysymbols, copy=False)
            else:
                expression = sp.sympify(expression)
                if expression.is_Pow:
                    assert len(expression.args) == 2
                    expression_base = Polynomial.from_expression(expression.args[0], polysymbols)
                    expression_exponent = Polynomial.from_expression(expression.args[1], polysymbols)
                    expression = ExponentiatedPolynomial(expression_base.expolist, expression_base.coeffs,
                                                         expression_exponent, polysymbols, copy=False)
                else:
                    expression = Polynomial.from_expression(expression, polysymbols)
                    expression = ExponentiatedPolynomial(expression.expolist, expression.coeffs,
                                                         Polynomial.from_expression(1, polysymbols), # exponent
                                                         polysymbols, copy=False)
        elif target_type == _Expression:
            expression = Expression(expression, polysymbols)
        else:
            raise RuntimeError('`_parse_expressions` is only implemented for `target_type`s in %s, not for %s. (raised while parsing `%s`)' \
                                % (set([_Expression, ExponentiatedPolynomial]), target_type, name_of_make_argument_being_parsed))
    parsed.append(expression)

def _convert_input(target_directory, name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   other_variables, form_optimization_level, form_work_space,
                   stabilize, contour_deformation_polynomial, decomposition_method):
    'Get the data types right.'

    # parse symbols
    integration_variables = sympify_symbols(list(integration_variables), 'All `integration_variables` must be symbols.')
    regulators = sympify_symbols(list(regulators), 'All `regulators` must be symbols.')
    polynomial_names = sympify_symbols(list(polynomial_names), 'All `polynomial_names` must be symbols.')

    # check format of `requested_orders` --> must be 1d and have the length as `regulators`
    requested_orders = np.array(requested_orders)
    assert len(requested_orders.shape) == 1, '`requested_orders` has the wrong shape. (is %s; should be (%i,))' % (requested_orders.shape,len(regulators))
    assert requested_orders.shape[0] == len(regulators), 'The length of `requested_orders` (%i) must match the length of `regulators` (%i)' % (len(requested_orders),len(regulators))


    # parse expressions
    # the `polynomials_to_decompose` should have ``integration_variables + regulators`` as polysymbols
    polynomials_to_decompose = _parse_expressions(polynomials_to_decompose, integration_variables + regulators, ExponentiatedPolynomial, 'polynomials_to_decompose')
    # the `other_polynomials` and the `remainder_expression` should additionally have the `polynomial_names` as polysymbols
    other_polynomials = _parse_expressions(other_polynomials, integration_variables + regulators + polynomial_names, ExponentiatedPolynomial, 'other_polynomials')
    remainder_expression = _parse_expressions(remainder_expression, integration_variables + regulators + polynomial_names, _Expression, 'remainder_expression')

    return (target_directory, name, integration_variables, regulators,
            requested_orders, polynomials_to_decompose, polynomial_names,
            other_polynomials, prefactor, remainder_expression, functions,
            other_variables, form_optimization_level, form_work_space,
            stabilize, contour_deformation_polynomial, decomposition_method)


# ---------------------------------- main function ----------------------------------
def make_package(target_directory, name, integration_variables, regulators, requested_orders,
                 polynomials_to_decompose, polynomial_names=[], other_polynomials=[],
                 prefactor=1, remainder_expression=1, functions=[], other_variables=[],
                 form_optimization_level=2, form_work_space='500M', stabilize=False,
                 contour_deformation_polynomial=None, decomposition_method='iterative_no_primary'):
    r'''
    Decompose, subtract and expand an expression.
    Return it as c++ package.

    :param target_directory:
        string;
        The output directory.

    :param name:
        string;
        The name of the c++ namepace.

    :param integration_variables:
        iterable of strings or sympy symbols;
        The variables that are to be integrated from
        0 to 1.

    :param regulators:
        iterable of strings or sympy symbols;
        The UV/IR regulators of the integral.

    :param requested_orders:
        iterable of integers;
        Compute the expansion in the regulators to these
        orders.

    :param polynomials_to_decompose:
        iterable of strings or sympy expressions or
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        or :class:`pySecDec.algebra.Polynomial`;
        The polynomials to be decomposed.

    :param polynomial_names:
        iterable of strings;
        Assign symbols for the `polynomials_to_decompose`.
        These can be referenced in the `other_polynomials`;
        see `other_polynomials` for details.

    :param other_polynomials:
        iterable of strings or sympy expressions or
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        or :class:`pySecDec.algebra.Polynomial`;
        Additional polynomials where no decomposition is
        attempted.
        The symbols defined in `polynomial_names`
        can be used to reference the
        `polynomials_to_decompose`. This is particularly
        useful when computing loop integrals where the
        "numerator" can depend on the first and second
        Symanzik polynomials.

        Example (1-loop bubble with numerator):

        >>> polynomials_to_decompose = ["(x0 + x1)**(2*eps - 4)",
        ...                             "(-p**2*x0*x1)**(-eps))"]
        >>> polynomial_names = ["U", "F"]
        >>> other_polynomials = ["""   (eps - 1)*s*U**2
        ...                          + (eps - 2)*F
        ...                          - (eps - 1)*2*s*x0*U
        ...                          + (eps - 1)*s*x0**2"""]

        .. seealso::
            :mod:`pySecDec.loop_integral`

    Note that the `polynomial_names` refer to the
    `polynomials_to_decompose` **without** their
    exponents.

    :param prefactor:
        string or sympy expression,
        optional;
        A factor that does not depend on the integration
        variables.

    :param remainder_expression:
        string or sympy expression or
        :class:`pySecDec.algebra._Expression`, optional;
        An additional factor.

        Dummy function must be provided with all arguments,
        e.g. ``remainder_expression='exp(eps)*f(x0,x1)'``.
        In addition, all dummy function must be listed
        in `functions`.

        The `polynomial_names` can be used like in
        `other_polynomials`. In particular, they go
        **without** arguments, e.g.
        ``remainder_expression='exp(g(U,F))'``.

    :param functions:
        iterable of strings or sympy symbols, optional;
        Function symbols occuring in `remainder_expression`,
        e.g.``[f]``.

    :param other_variables:
        iterable of strings or sympy symbols, optional;
        Other occuring symbols.

    :param form_optimization_level:
        integer out of the interval [0,3], optional;
        The optimization level to be used in FORM.
        Default: ``2``.

    :param form_work_space:
        string, optional;
        The FORM WorkSpace. Default: ``'500M'``.

    :param stabilize:
        bool, optional;
        Whether or not to bring subexpression to a common
        denominator. Default: ``False``.

        .. warning::
            This is very extensive concerning both - the
            algebra and the numerics. It should only be
            set to ``True`` if numerical instabilities
            occur.

    :param contour_deformation_polynomial:
        string or sympy symbol, optional;
        The name of the polynomial in
        `polynomials_to_decompose` or `other_polynomials`
        that is to be continued to the complex plane
        according to a :math:`- i\delta` prescription.
        For loop integrals, this is the second Symanzik
        polynomial ``F``.
        If not provided, no code for contour deformation
        is created.

    :param decomposition_method:
        string, optional;
        The strategy to decompose the polynomials. The
        following strategies are available:

        * 'iterative_no_primary' (default)
        * 'iterative'
        * 'geometric'

        'iterative' and 'geometric' are only valid for
        loop integrals. An end user should always use
        the default 'iterative_no_primary'.

    '''
    # convert input data types to the data types we need
    target_directory, name, integration_variables, regulators,
    requested_orders, polynomials_to_decompose, polynomial_names,
    other_polynomials, prefactor, remainder_expression, functions,
    other_variables, form_optimization_level, form_work_space,
    stabilize, contour_deformation_polynomial, decomposition_method = \
    _convert_input(target_directory, name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   other_variables, form_optimization_level, form_work_space,
                   stabilize, contour_deformation_polynomial, decomposition_method)

    # TODO: implement actual computation and code writing
