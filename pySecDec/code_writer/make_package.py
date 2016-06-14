"""
This module implements the main program - the function
:func:`.make_package`.

"""

from ..misc import sympify_symbols
from ..algebra import _Expression, Expression, Polynomial, \
                      ExponentiatedPolynomial, Pow, Product, \
                      DerivativeTracker, Function, Sum
from .. import decomposition
from ..subtraction import integrate_pole_part
from ..expansion import expand_singular, expand_Taylor
from ..misc import lowest_order
from .template_parser import parse_template_file, parse_template_tree
import numpy as np
import sympy as sp
from itertools import chain
import os

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
    return parsed

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
    other_variables = sympify_symbols(list(other_variables), 'All `other_variables` must be symbols.')
    functions = sympify_symbols(list(functions), 'All `functions` must be symbols.')
    if contour_deformation_polynomial is not None:
        contour_deformation_polynomial = sympify_symbols([contour_deformation_polynomial], '`contour_deformation_polynomial` must be a symbol.')[0]

    # define the symbols of the different classes of `_Expression`s
    symbols_polynomials_to_decompose = integration_variables + regulators
    all_symbols = symbols_other_polynomials = symbols_remainder_expression = integration_variables + regulators + polynomial_names

    # check format of `requested_orders` --> must be 1d and have the length as `regulators`
    requested_orders = np.array(requested_orders)
    assert len(requested_orders.shape) == 1, '`requested_orders` has the wrong shape. (is %s; should be (%i,))' % (requested_orders.shape,len(regulators))
    assert requested_orders.shape[0] == len(regulators), 'The length of `requested_orders` (%i) must match the length of `regulators` (%i)' % (len(requested_orders),len(regulators))

    # parse expressions
    polynomials_to_decompose = _parse_expressions(polynomials_to_decompose, symbols_polynomials_to_decompose, ExponentiatedPolynomial, 'polynomials_to_decompose')
    other_polynomials = _parse_expressions(other_polynomials, symbols_other_polynomials, ExponentiatedPolynomial, 'other_polynomials')
    remainder_expression = _parse_expressions([remainder_expression], symbols_remainder_expression, _Expression, 'remainder_expression')[0]

    # convert ``prefactor`` to sympy expression
    prefactor = sp.sympify(prefactor)

    # convert ``requested_orders`` to numpy array
    requested_orders = np.array(requested_orders)

    return (target_directory, name, integration_variables, regulators,
            requested_orders, polynomials_to_decompose, polynomial_names,
            other_polynomials, prefactor, remainder_expression, functions,
            other_variables, form_optimization_level, form_work_space,
            stabilize, contour_deformation_polynomial, decomposition_method,
            symbols_polynomials_to_decompose, symbols_other_polynomials,
            symbols_remainder_expression, all_symbols)


# ---------------------------------- decomposition ----------------------------------
_decomposition_strategies = dict(
                                    iterative=           dict(
                                                                 primary=decomposition.iterative.primary_decomposition,
                                                                 secondary=decomposition.iterative.iterative_decomposition
                                                             ),
                                    geometric=           dict(
                                                                 primary=decomposition.iterative.primary_decomposition,
                                                                 secondary=decomposition.iterative.iterative_decomposition
                                                             ),
                                    iterative_no_primary=dict(
                                                                 primary=lambda x: [x], # no primary decomposition
                                                                 secondary=decomposition.iterative.iterative_decomposition
                                                             )
                                )


# -------------------------------- template parsing ---------------------------------
def _parse_global_templates(target_directory, name, integration_variables,
                            regulators, polynomial_names, functions,
                            other_variables, form_optimization_level,
                            form_work_space, stabilize, requested_orders,
                            contour_deformation_polynomial):
    '''
    Create the `target_directory` and return the two
    optional arguments passed to :func:`parse_template_tree`.

    '''
    # TODO: check validity of symbol names and `name` in FORM and c++ (FORM: no underscores; both: no special characters, no leading numbers), and its length (due to FORM's line breaks)

    # initialize template replacements
    template_replacements = dict(
                                     name = name,
                                     integration_variables = _make_FORM_list(integration_variables),
                                     regulators = _make_FORM_list(regulators),
                                     polynomial_names = _make_FORM_list(regulators),
                                     functions = _make_FORM_list(functions),
                                     other_variables = _make_FORM_list(other_variables),
                                     form_optimization_level = form_optimization_level,
                                     form_work_space = form_work_space,
                                     contour_deformation = int(contour_deformation_polynomial is not None),
                                     stabilize = int(stabilize),
                                     requested_orders = _make_FORM_list(requested_orders)
                                )

    # configure template parser
    file_renamings = {
                          # replace "name" by the name of the integral
                          'name' : name,
                          'name.hpp' : name + '.hpp',

                          # the files below are specific for each sector --> do not parse globally
                          'contour_deformation.h' : None,
                          'sector.h' : None,

                          # "integrands.hpp" can only be written after the decomposition is completed
                          'integrands.hpp' : None
                     }

    # the files below are only relevant for contour deformation --> do not parse if deactivated
    if contour_deformation_polynomial is not None:
        for filename in ['write_contourdef.frm', 'optimize.hpp', 'contour_deformation.h']:
            file_renamings[filename] = None

    # get path to the directory with the template files (path relative to directory with this file: "./templates/")
    from . import test_make_package as _unittests_for_this_module
    template_sources = os.path.join(os.path.split(_unittests_for_this_module.__file__)[0],'templates')

    # initialize the target directory with the sector independent files
    parse_template_tree(template_sources, target_directory, template_replacements, file_renamings)

    # return parser options
    return template_sources, target_directory, template_replacements, file_renamings


# --------------------------------- write FORM code ---------------------------------
def _make_FORM_list(python_list):
    '''
    Convert a python list to a string to be used
    in FORM like a list.

    Example: ``['a', 'b', 'c'] --> 'a, b, c'``

    '''
    return ', '.join(str(item) for item in python_list)

def _make_FORM_definition(name, expression):
    r'''
    Write the following line for the insertion
    in FORM:

    #define `name` "`expression`"

    Sample output:
    #define f " + (3)*y*z + (1)*x^2"

    '''
    return '#define %s "%s"\n' % \
    (
        name,
        str(expression).replace('**','^')
    )

def _make_FORM_shifted_orders(positive_powers):
    r'''
    Write FORM code that defines the preprocessor
    variables
    `shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex''.

    '''
    codelines = []
    orders_index = 0
    for multiindex in positive_powers:
        orders_index += 1
        regulator_index = 0
        for regulator_power in multiindex:
            regulator_index += 1
            codelines.append( '#define shiftedRegulator%iPowerOrder%i "%i"' % (regulator_index,orders_index,regulator_power) )
    return '\n'.join(codelines)

def _derivative_muliindex_to_name(basename, multiindex):
    '''
    Convert a derivative multiindex as returned by
    :meth:`pySecDec.algebra.DerivativeTracker.compute_derivatives`
    to the form ``d...d<basename>d<index>...d<index>``.

    Example:

    >>> _derivative_muliindex_to_name('f', (1,2,1))
    ddddfd0d1d1d2

    '''
    prefix = 'd' * sum(multiindex)
    suffix = ''.join(('d' + str(i)) * depth for i,depth in enumerate(multiindex))
    return prefix + basename + suffix

# define the internal names to be used in FORM
internal_prefix = 'SecDecInternal'
FORM_names = dict(
    cal_I=internal_prefix+'CalI',
    regular=internal_prefix+'Regular'
)


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
        e.g.``['f']``.

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
        The name of the polynomial in `polynomial_names`
        that is to be continued to the complex plane
        according to a :math:`- i\delta` prescription.
        For loop integrals, this is the second Symanzik
        polynomial ``F``.
        If not provided, no code for contour deformation
        is created.

        .. note::
            The `contour_deformation_polynomial` must **NOT**
            depend on the `regulators`.

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
    target_directory, name, integration_variables, regulators, \
    requested_orders, polynomials_to_decompose, polynomial_names, \
    other_polynomials, prefactor, remainder_expression, functions, \
    other_variables, form_optimization_level, form_work_space, \
    stabilize, contour_deformation_polynomial, decomposition_method, \
    symbols_polynomials_to_decompose, symbols_other_polynomials, \
    symbols_remainder_expression, all_symbols = \
    _convert_input(target_directory, name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   other_variables, form_optimization_level, form_work_space,
                   stabilize, contour_deformation_polynomial, decomposition_method)

    # configure the template parser and parse global files
    template_sources, target_directory, template_replacements, file_renamings = \
        _parse_global_templates(
        target_directory, name, integration_variables, regulators,
        polynomial_names, functions, other_variables,
        form_optimization_level, form_work_space, stabilize,
        requested_orders, contour_deformation_polynomial
    )

    # get the highest poles from the ``prefactor``
    highest_prefactor_pole_orders = np.array([lowest_order(prefactor, regulator) for regulator in regulators])

    # compute the required expansion order accounting for the prefactor
    required_orders = requested_orders + highest_prefactor_pole_orders

    # get the decomposition routines
    strategy = _decomposition_strategies[decomposition_method]

    # define the `Polynomial` "x0 + x1 + x2 + ..." to keep track of the transformations
    transformations = Polynomial(np.identity(len(integration_variables), dtype=int), [1]*len(integration_variables), integration_variables)

    # hide ``regulators`` and ``polynomial_names`` from decomposition
    polynomials_to_decompose_hidden_regulators = []
    polynomials_to_decompose_hide_containers = []
    for poly in polynomials_to_decompose:
        poly, hidden = decomposition.hide(poly, len(regulators))
        polynomials_to_decompose_hidden_regulators.append(poly)
        polynomials_to_decompose_hide_containers.append(hidden)
    other_polynomials_hidden_regulators_and_names = []
    other_polynomials_name_hide_containers = []
    other_polynomials_regulator_hide_containers = []
    for poly in other_polynomials:
        poly, hidden = decomposition.hide(poly, len(polynomial_names))
        other_polynomials_name_hide_containers.append(hidden)
        poly, hidden = decomposition.hide(poly, len(regulators))
        other_polynomials_regulator_hide_containers.append(hidden)
        other_polynomials_hidden_regulators_and_names.append(poly)

    # initialize the decomposition
    initial_sector = decomposition.Sector(polynomials_to_decompose_hidden_regulators, other_polynomials_hidden_regulators_and_names + [transformations])

    # initialize the counter
    sector_index = 0

    # define the imaginary unit
    imaginary_unit = sp.sympify('I') # TODO: replace by "i_" in FORM

    # we must backwards traverse the `polynomial_names` --> create the reversed list once and for all
    reversed_polynomial_names = list(polynomial_names) # copy
    reversed_polynomial_names.reverse()

    # get the index of the `contour_deformation_polynomial` in `polynomial_names`
    if contour_deformation_polynomial is not None:
        str_contour_deformation_polynomial = str(contour_deformation_polynomial)
        contour_deformation_polynomial_index = 0
        while str(polynomial_names[contour_deformation_polynomial_index]) != str_contour_deformation_polynomial:
            contour_deformation_polynomial_index += 1
            if len(polynomial_names) <= contour_deformation_polynomial_index:
                raise IndexError('Could not find the `contour_deformation_polynomial` "%s" in `polynomial_names`.' % str_contour_deformation_polynomial)

    for primary_sector in strategy['primary'](initial_sector):

        # primary decomposition removes one integration parameter --> redefine `integration_variables` and the symbols of the different classes of `_Expression`s
        integration_variables = list(primary_sector.Jacobian.polysymbols) # make a copy
        symbols_polynomials_to_decompose = integration_variables + regulators
        all_symbols = symbols_other_polynomials = symbols_remainder_expression = integration_variables + regulators + polynomial_names

        # get the indices of the `integration_variables`, the `regulators`, and the `polynomial_names` in the polysymbols
        integration_variable_indices = list(range(len(integration_variables)))
        regulator_indices = [i + len(integration_variables) for i in range(len(regulators))]

        # define "elementary" `_Expression`s such as ``x = Polynomial.from_expression('x', polysymbols)`` for all x
        elementary_monomials = []
        for i, symbol in enumerate(symbols_polynomials_to_decompose):
            expolist = np.zeros([1,len(symbols_polynomials_to_decompose)], dtype=int)
            expolist[:,i] = 1
            elementary_monomials.append( Polynomial(expolist, np.array([1]), symbols_polynomials_to_decompose, copy=False) )

        if contour_deformation_polynomial is not None: # TODO: further checks on `contour_deformation_polynomial`, e.g. it must not depend on the regulators
            # Need all first and second derivatives of the `contour_deformation_polynomial`.
            # Since the `contour_deformation_polynomial` is left symbolic they are equal for every subsector after primary decomposition.
            symbolic_contour_deformation_polynomial = Function(str(contour_deformation_polynomial), *symbols_polynomials_to_decompose) # TODO: check that the referenced polynomial does not depend on the `regulators`
            symbolic_contour_deformation_polynomial = DerivativeTracker(symbolic_contour_deformation_polynomial)

            # compute the transformation of the integration parameters and its Jacobian matrix (see e.g. section 3.2 in arXiv:1601.03982):
            # ``z_k({x_k}) = x_k - i * lambda_k * (1 - x_k) * Re(dF_dx_k)``, where "dF_dx_k" denotes the derivative of ``F`` by ``x_k``
            # Remarks:
            #   - We account for the Re(...) in the numerics only.
            #   - The determinant of the Jacobian matrix is calculated numerically.
            deformation_parameters = [sp.symbols('lambda%i'%i) for i in range(len(integration_variables))] # TODO: make these symbols user defined
            transformed_integration_parameters = [
                                                     integration_variables[k] +
                                                     ( - imaginary_unit * deformation_parameters[k] * integration_variables[k] * (1 - integration_variables[k]) ) *
                                                     symbolic_contour_deformation_polynomial.simplify().derive(k)
                                                     for k in range(len(integration_variables))
                                                 ]

            contourdef_Jacobian = np.empty((len(integration_variables), len(integration_variables)), dtype=object)
            for i in range(len(integration_variables)):
                for j in range(len(integration_variables)):
                    contourdef_Jacobian[i,j] = transformed_integration_parameters[i].simplify().derive(j)

        # insert the `polynomials_to_decompose` as dummy functions into `other_polynomials` and `remainder_expression`
        # we want to remove them from the symbols --> traverse backwards and pop the last of the `polysymbols`
        # redefine ``all_symbols`` (we are going to removed some)
        all_symbols = integration_variables + regulators
        this_sector_remainder_expression = remainder_expression
        for i in range(len(primary_sector.other) - 1): # "-1" because of ``transformations``
            decomposition.unhide(primary_sector.other[i], other_polynomials_name_hide_containers[i])
            for poly_name in reversed_polynomial_names:
                primary_sector.other[i] = primary_sector.other[i].replace(-1, poly_name(*all_symbols), remove=True)
                this_sector_remainder_expression = this_sector_remainder_expression.replace(-1, poly_name(*all_symbols), remove=True)

        # convert the coefficients and exponents to `_Expression`s
        # TODO: implement this step, the polysymbols should be `all_symbols`

        # we later need ``1`` packed into specific types
        polynomial_one = Polynomial(np.zeros([1,len(all_symbols)], dtype=int), np.array([1]), all_symbols, copy=False)
        pole_part_initializer = Pow(polynomial_one, -polynomial_one)

        # TODO: exponents should not depend on the `integration_variables` --> assert in `_convert_input`

        for sector in strategy['secondary'](primary_sector):
            sector_index += 1

            # unhide the regulator
            #  - in the Jacobian
            Jacobian = sector.Jacobian
            Jacobian.number_of_variables += len(regulators)
            Jacobian.polysymbols += regulators
            Jacobian.expolist = np.hstack([Jacobian.expolist, np.zeros([len(Jacobian.coeffs),len(regulators)], dtype=int)])

            #  - in ``sector.cast``
            for product,hidden_regulators in zip(sector.cast,polynomials_to_decompose_hide_containers):
                mono, poly = product.factors
                mono.number_of_variables += len(regulators)
                mono.polysymbols += regulators
                mono.expolist = np.hstack([mono.expolist, np.zeros([1,len(regulators)], dtype=int)])
                decomposition.unhide(poly, hidden_regulators)

            #  - in ``sector.other``
            for i in range(len(sector.other) - 1): # "-1" because of ``transformations``
                decomposition.unhide(sector.other[i], other_polynomials_regulator_hide_containers[i])

            #  - in the ``transformations``
            this_transformation = sector.other.pop() # remove transformation from ``sector.other``
            this_transformation.number_of_variables += len(regulators)
            this_transformation.polysymbols += regulators
            this_transformation.expolist = np.hstack([this_transformation.expolist, np.zeros([len(this_transformation.coeffs),len(regulators)], dtype=int)])


            # insert the monomials of the `polynomial_names` into ``sector.other`` ``<symbol> --> xi**pow_i * <symbol>``
            # ``[:,:-len(regulators)]`` is the part of the expolist that contains only the powers of the `integration_variables` but not of the `regulators`;
            #     i.e. it excludes the powers of the regulators from the modification
            for poly in sector.other:
                for name,(mono,_),hidden_name in zip(polynomial_names,sector.cast,other_polynomials_name_hide_containers):
                    poly.expolist[:,:-len(regulators)] += np.einsum('i,k->ik', hidden_name.expolist[:,0], mono.expolist[0,:-len(regulators)])

            # factorize polynomials in ``sector.other``
            for i,to_factorize in enumerate(sector.other):
                to_factorize = sector.other[i] = Product(polynomial_one, to_factorize)
                decomposition.refactorize(to_factorize)

            # subtraction needs type `ExponentiatedPolynomial` for all factors in its monomial part
            Jacobian = ExponentiatedPolynomial(Jacobian.expolist, Jacobian.coeffs, polysymbols=Jacobian.polysymbols, exponent=polynomial_one, copy=False)

            # initialize the product of monomials for the subtraction
            monomial_factors = chain([Jacobian], (prod.factors[0] for prod in sector.cast), (prod.factors[0] for prod in sector.other))
            monomials = Product(*monomial_factors, copy=False)

            # define ``cal_I``, the part of the integrand that does not lead to poles
            # TODO: apply ``transformations`` to ``this_sector_remainder_expression`` BEFORE taking derivatives
            cal_I = Product(this_sector_remainder_expression, *(prod.factors[1] for prod in chain(sector.cast, sector.other)), copy=False)

            # it is faster to use a dummy function for ``cal_I`` and substitute back in FORM
            # symbolic cal_I wrapped in a `DerivativeTracker` to keep track of derivatives
            symbolic_cal_I = DerivativeTracker(Function(FORM_names['cal_I'], *elementary_monomials))

            # initialize the Product to be passed to `integrate_pole_part` (the subtraction) and subtract
            subtraction_initializer = Product(monomials, pole_part_initializer, symbolic_cal_I, copy=False)
            subtracted = integrate_pole_part(subtraction_initializer, *integration_variable_indices)

            # intialize expansion
            pole_parts = [s.factors[1].simplify() for s in subtracted]
            regular_parts = [Product( *([s.factors[0]] + s.factors[2:]), copy=False ) for s in subtracted]
            # introduce dummy functions --> faster in python
            symbolic_regular_parts = \
            [
                DerivativeTracker\
                (
                    Function\
                    (
                        FORM_names['regular'] + str(i), *elementary_monomials
                    ), copy=False
                )
                for i,_ in enumerate(subtracted)
            ]

            # expand poles
            integrand_summands = []
            for i,(regular,singular) in enumerate(zip(symbolic_regular_parts, pole_parts)):
                # must expand every term to the requested order plus the highest pole it multiplies
                # We calculated the highest pole order of the prefactor (variable ``highest_prefactor_pole_orders``) above.
                # In addition, we have to take the poles of the current term into account.
                singular_expanded = expand_singular(Product(singular, copy=False), regulator_indices, required_orders)

                highest_poles_current_term = - singular_expanded.expolist[:,regulator_indices].min(axis=0)
                expansion_orders = requested_orders + highest_prefactor_pole_orders + highest_poles_current_term
                regular_expanded = expand_Taylor(regular, regulator_indices, expansion_orders)

                if i == 0: # first iteration; ``highest_poles_current_sector`` not yet set
                    highest_poles_current_sector = highest_poles_current_term
                else:
                    highest_poles_current_sector = np.maximum(highest_poles_current_sector, highest_poles_current_term)

                # TODO: the following multiplication can generate terms of higher orders than requested --> discard them here or in FORM?
                integrand_summands.append(singular_expanded * regular_expanded)

            integrand = Sum(*integrand_summands, copy=False)

            # compute the required derivatives
            derivatives = {}
            def update_derivatives(basename, derivative_tracker, full_expression):
                derivatives[basename] = full_expression # include undifferentiated expression
                for multiindex,expression in derivative_tracker.compute_derivatives(full_expression).items():
                    derivatives[_derivative_muliindex_to_name(basename, multiindex)] = expression

            #  - for the contour deformation
            if contour_deformation_polynomial is not None:
                update_derivatives(
                    contour_deformation_polynomial, # basename
                    symbolic_contour_deformation_polynomial, # derivative tracker
                    sector.cast[contour_deformation_polynomial_index] # full expression
                )

            #  - for the "regular parts" arising in the subtraction
            for index, (regular_part, symbolic_regular_part) in enumerate(zip(regular_parts, symbolic_regular_parts)):
                update_derivatives(basename=FORM_names['regular']+str(index), derivative_tracker=symbolic_regular_part, full_expression=regular_part)

            #  - for the part of the integrand that does not lead to poles
            # Must inherit derivatives from the "regular parts" because the derivative
            # tracker does not work any more after "cal_I" is hidden in them.
                symbolic_cal_I.derivatives.update(symbolic_regular_part.derivatives)
            update_derivatives(basename=FORM_names['cal_I'], derivative_tracker=symbolic_cal_I, full_expression=cal_I)


            # generate the function definitions the insertion in FORM
            form_function_definitions = ''.join(
                _make_FORM_definition(name, expression)
                for name, expression in derivatives.items()
            )
            form_insertions = _make_FORM_list(derivatives.keys())

            print(form_function_definitions)
            print(form_insertions)
            #for i,p in enumerate(symbolic_regular_parts):
            #    print(p.compute_derivatives(regular_parts[i]))
            #    print()
            #    print()

            #print(integrand)

# TODO: compute required derivatives
# TODO: implement code writing
