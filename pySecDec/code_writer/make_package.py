"""
This module implements the main program - the function
:func:`.make_package`.

"""

from __future__ import print_function
from ..misc import sympify_symbols, rangecomb
from ..algebra import _Expression, Expression, Polynomial, \
                      ExponentiatedPolynomial, Pow, Product, \
                      ProductRule, DerivativeTracker, Function, \
                      Sum
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

_sympy_one = sp.sympify(1)

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
                expression.exponent = sp.sympify(expression.exponent)
                expression.coeffs = np.array( sp.sympify(coeff) for coeff in expression.coeffs )
            elif type(expression) == Polynomial:
                expression = ExponentiatedPolynomial(expression.expolist,
                                                     np.array( sp.sympify(coeff) for coeff in expression.coeffs ),
                                                     _sympy_one, # exponent
                                                     polysymbols, copy=False)
            else:
                expression = sp.sympify(expression)
                if expression.is_Pow:
                    assert len(expression.args) == 2
                    expression_base = Polynomial.from_expression(expression.args[0], polysymbols)
                    expression_exponent = sp.sympify(expression.args[1], polysymbols)
                    expression = ExponentiatedPolynomial(expression_base.expolist,
                                                         expression_base.coeffs,
                                                         expression_exponent, polysymbols, copy=False)
                else:
                    expression = Polynomial.from_expression(expression, polysymbols)
                    expression = ExponentiatedPolynomial(expression.expolist,
                                                         expression.coeffs,
                                                         _sympy_one, # exponent
                                                         polysymbols, copy=False)
        elif target_type == _Expression:
            if not isinstance(expression, _Expression):
                expression = Expression(expression, polysymbols)
        else:
            raise RuntimeError('`_parse_expressions` is only implemented for `target_type`s in %s, not for %s. (raised while parsing `%s`)' \
                                % (set([_Expression, ExponentiatedPolynomial]), target_type, name_of_make_argument_being_parsed))
        parsed.append(expression)
    return parsed

def _convert_input(target_directory, name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   real_parameters, complex_parameters, form_optimization_level,
                   form_work_space, form_insertion_depth, stabilize,
                   contour_deformation_polynomial, decomposition_method):
    'Get the data types right.'

    # parse symbols
    integration_variables = sympify_symbols(list(integration_variables), 'All `integration_variables` must be symbols.')
    regulators = sympify_symbols(list(regulators), 'All `regulators` must be symbols.')
    polynomial_names = sympify_symbols(list(polynomial_names), 'All `polynomial_names` must be symbols.')
    real_parameters= sympify_symbols(list(real_parameters), 'All `real_parameters` must be symbols.')
    complex_parameters = sympify_symbols(list(complex_parameters), 'All `complex_parameters` must be symbols.')
    functions = sympify_symbols(list(functions), 'All `functions` must be symbols.')
    if contour_deformation_polynomial is not None:
        contour_deformation_polynomial = sympify_symbols([contour_deformation_polynomial], '`contour_deformation_polynomial` must be a symbol.')[0]

    # define the symbols of the different classes of `_Expression`s
    symbols_polynomials_to_decompose = integration_variables + regulators
    symbols_remainder_expression = integration_variables + polynomial_names
    all_symbols = symbols_other_polynomials = integration_variables + regulators + polynomial_names

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

    # convert ``form_insertion_depth`` to integer and check nonnegativity
    assert form_insertion_depth == int(form_insertion_depth), '`form_insertion_depth` must be an integer.'
    assert form_insertion_depth >= 0, '`form_insertion_depth` must not be negative.'

    return (target_directory, name, integration_variables, regulators,
            requested_orders, polynomials_to_decompose, polynomial_names,
            other_polynomials, prefactor, remainder_expression, functions,
            real_parameters, complex_parameters, form_optimization_level,
            form_work_space, form_insertion_depth, stabilize,
            contour_deformation_polynomial, decomposition_method,
            symbols_polynomials_to_decompose, symbols_other_polynomials,
            symbols_remainder_expression, all_symbols)


# ---------------------------------- decomposition ----------------------------------
_decomposition_strategies = dict(
                                    iterative=           dict(
                                                                 primary=decomposition.iterative.primary_decomposition,
                                                                 secondary=decomposition.iterative.iterative_decomposition
                                                             ),
                                    geometric=           dict(
                                                                 primary=lambda x: [decomposition.geometric.Cheng_Wu(x)],
                                                                 secondary=decomposition.geometric.geometric_decomposition
                                                             ),
                                    iterative_no_primary=dict(
                                                                 primary=lambda x: [x], # no primary decomposition
                                                                 secondary=decomposition.iterative.iterative_decomposition
                                                             )
                                )


# -------------------------------- template parsing ---------------------------------
def _parse_global_templates(target_directory, name, regulators, polynomial_names,
                            real_parameters, complex_parameters, form_optimization_level,
                            form_work_space, form_insertion_depth, stabilize,
                            requested_orders, contour_deformation_polynomial,
                            integrand_container_type):
    '''
    Create the `target_directory` and return the two
    optional arguments passed to :func:`parse_template_tree`.

    '''
    # TODO: check validity of symbol names and `name` in FORM and c++ (FORM: no underscores; both: no special characters, no leading numbers)

    # initialize template replacements
    template_replacements = dict(
                                     name = name,
                                     number_of_real_parameters = len(real_parameters),
                                     real_parameters = _make_FORM_list(real_parameters),
                                     number_of_complex_parameters = len(complex_parameters),
                                     complex_parameters = _make_FORM_list(complex_parameters),
                                     have_complex_parameters = len(complex_parameters) > 0,
                                     number_of_regulators = len(regulators),
                                     regulators = _make_FORM_list(regulators),
                                     polynomial_names = _make_FORM_list(regulators),
                                     form_optimization_level = form_optimization_level,
                                     form_work_space = form_work_space,
                                     form_insertion_depth = form_insertion_depth,
                                     contour_deformation = int(contour_deformation_polynomial is not None),
                                     stabilize = int(stabilize),
                                     requested_orders = _make_FORM_list(requested_orders),
                                     integrand_container_type = integrand_container_type
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
    if contour_deformation_polynomial is None:
        for filename in ['write_contourdef.frm', 'contour_deformation.h']:
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

    Example: ``['a', 'b', 'c'] --> 'a,b,c'``

    '''
    return ','.join(str(item) for item in python_list)

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

def _make_FORM_function_definition(name, expression, limit):
    '''
    Split an `expression` whose string representation
    exceeds the `limit` of characters. Return a FORM
    procedure "defineReplacementExpression`name`"
    that defines `expression` as local expression
    named ``replacement``.

    :param name:
        string;
        The name of the expression. It is part of the
        name of the resulting FORM procedure
        "defineReplacementExpression`name`".

    :param expression:
        :class:`pySecDec.algebra._Expression`;
        The expression to be split.

    :param limit:
        integer;
        The maximum number of characters in each
        subexpression.

    '''
    first_line = "#procedure generateReplacement%s(?replaceArg)" % name
    last_line = "#endProcedure"
    codelines = []

    def recursion(start_i, expression):
        str_expression = str(expression)

        if len(str_expression) <= limit: # no need to split
            codelines.append( "  L tmp%i = (%s)*replace_(`?replaceArg');" % (start_i,str_expression) )
            return start_i + 1

        if type(expression) == ProductRule:
            expression = expression.to_sum()

        if type(expression) == Sum:
            relevant_indices = []
            i = start_i
            for summand in expression.summands:
                i = recursion(i, summand)
                relevant_indices.append(i - 1)
            codelines.append( "  .sort" )
            codelines.append( "  L tmp%i = " % i + "+".join("tmp%i" % j for j in relevant_indices) + ";" )
            return i + 1

        if type(expression) == Product:
            relevant_indices = []
            i = start_i
            for factor in expression.factors:
                i = recursion(i, factor)
                relevant_indices.append(i - 1)
            codelines.append( "  .sort" )
            codelines.append( "  L tmp%i = " % i + "*".join("tmp%i" % j for j in relevant_indices) + ";" )
            return i + 1

        # rescue: print warning and write unsplit expression
        print( 'WARNING: Could not split "%s" (not implemented for %s)' % (name,type(expression)) )
        codelines.append("  L tmp%i = (%s)*replace_(`?replaceArg');" % (start_i, str_expression))
        return start_i + 1

    largest_i = recursion(1, expression) - 1
    codelines.append( "  .sort" )
    codelines.append( "  drop " + ",".join("tmp%i" % j for j in range(1,largest_i+1)) + ";" )
    codelines.append( "  L replacement = tmp%i;" % largest_i )
    codelines.append( "  .sort" )
    return "\n".join(chain([first_line],codelines,[last_line,'']))

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

def _make_FORM_Series_initilization(min_orders, max_orders, sector_ID):
    '''
    Write the c++ code that initilizes the container class
    (``Series<Series<...<Series<IntegrandContainer>>...>``).

    '''
    assert len(min_orders) == len(max_orders)
    last_regulator_index = len(min_orders) - 1

    def multiindex_to_cpp_order(multiindex):
        '(-1,3,2,-4) --> n1_3_2_n4'
        snippets = []
        for order in multiindex:
            snippets.append(str(order).replace('-','n'))
        return '_'.join(snippets)

    current_orders = np.array(min_orders) # use as nonlocal variable in `recursion`
    def recursion(regulator_index):
        if regulator_index < last_regulator_index:
            outstr_body_snippets = []
            outstr_head = '{%i,%i,{' % (min_orders[regulator_index],max_orders[regulator_index])
            for this_regulator_order in range(min_orders[regulator_index],max_orders[regulator_index]+1):
                current_orders[regulator_index] = this_regulator_order
                outstr_body_snippets.append( recursion(regulator_index + 1) )
            outstr_tail = '},true}'
            return ''.join( (outstr_head, ','.join(outstr_body_snippets), outstr_tail) )
        else: # regulator_index == last_regulator_index; i.e. processing last regulator
            outstr_head = '{%i,%i,{{' % (min_orders[regulator_index],max_orders[regulator_index])
            outstr_body_snippets = []
            for this_regulator_order in range(min_orders[regulator_index],max_orders[regulator_index]+1):
                current_orders[regulator_index] = this_regulator_order
                outstr_body_snippets.append(
                    'sector_%(sector_ID)i_order_%(cpp_order)s_numIV,sector_%(sector_ID)i_order_%(cpp_order)s_integrand' \
                    % dict(sector_ID=sector_ID,cpp_order=multiindex_to_cpp_order(current_orders))
                )
            outstr_tail = '}},true}'
            return ''.join( (outstr_head, '},{'.join(outstr_body_snippets), outstr_tail) )

    return recursion(0)


# ---------------------------------- main function ----------------------------------
def make_package(target_directory, name, integration_variables, regulators, requested_orders,
                 polynomials_to_decompose, polynomial_names=[], other_polynomials=[],
                 prefactor=1, remainder_expression=1, functions=[], real_parameters=[],
                 complex_parameters=[], form_optimization_level=2, form_work_space='500M',
                 form_insertion_depth=1, stabilize=False, contour_deformation_polynomial=None,
                 decomposition_method='iterative_no_primary'):
    r'''
    Decompose, subtract and expand an expression.
    Return it as c++ package.

    .. note::
        Use ``I`` to denote the imaginary unit.

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

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as real variables.

    :param complex_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as complex variables.

    :param form_optimization_level:
        integer out of the interval [0,3], optional;
        The optimization level to be used in FORM.
        Default: ``2``.

    :param form_work_space:
        string, optional;
        The FORM WorkSpace. Default: ``'500M'``.

    :param form_insertion_depth:
        nonnegative integer, optional;
        How deep FORM should try resolving nested function
        calls. Default: ``1``.

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
    real_parameters, complex_parameters, form_optimization_level, \
    form_work_space, form_insertion_depth, stabilize, \
    contour_deformation_polynomial, decomposition_method, \
    symbols_polynomials_to_decompose, symbols_other_polynomials, \
    symbols_remainder_expression, all_symbols = \
    _convert_input(target_directory, name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   real_parameters, complex_parameters, form_optimization_level,
                   form_work_space, form_insertion_depth, stabilize,
                   contour_deformation_polynomial, decomposition_method)

    # construct the c++ type of the integrand container class
    # for two regulators, the resulting code should read:
    # "secdecutil::Series<secdecutil::Series<IntegrandContainer>>"
    integrand_container_type = 'secdecutil::Series<' * len(regulators) + 'IntegrandContainer' + '>' * len(regulators)

    # configure the template parser and parse global files
    template_sources, target_directory, template_replacements, file_renamings = \
        _parse_global_templates(
        target_directory, name, regulators, polynomial_names,
        real_parameters, complex_parameters, form_optimization_level,
        form_work_space, form_insertion_depth, stabilize,
        requested_orders, contour_deformation_polynomial,
        integrand_container_type
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
    imaginary_unit = sp.sympify('I')

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
        symbols_polynomials_to_decompose = symbols_other_polynomials = integration_variables + regulators
        symbols_remainder_expression = integration_variables + polynomial_names
        all_symbols = integration_variables + regulators + polynomial_names

        # define `integration_variables` in the template system
        template_replacements['number_of_integration_variables'] = len(integration_variables)
        template_replacements['integration_variables'] = _make_FORM_list(integration_variables)

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
        # redefine ``all_symbols`` and ``symbols_remainder_expression``
        all_symbols = symbols_remainder_expression = integration_variables + regulators
        this_primary_sector_remainder_expression = remainder_expression
        for i in range(len(primary_sector.other) - 1): # "-1" because of ``transformations``
            poly = primary_sector.other[i]
            decomposition.unhide(poly, other_polynomials_name_hide_containers[i])
            for poly_name in reversed_polynomial_names:
                poly = poly.replace(-1, poly_name(*all_symbols), remove=True)
            primary_sector.other[i] = poly
        for poly_name in reversed_polynomial_names:
            this_primary_sector_remainder_expression = this_primary_sector_remainder_expression.replace(-1, poly_name(*all_symbols), remove=True)

        # we later need ``1`` packed into specific types
        polynomial_one = Polynomial(np.zeros([1,len(all_symbols)], dtype=int), np.array([1]), all_symbols, copy=False)
        pole_part_initializer = Pow(polynomial_one, -polynomial_one)

        # TODO: exponents should not depend on the `integration_variables` --> assert in `_convert_input`

        for sector in strategy['secondary'](primary_sector):
            sector_index += 1

            # extract ``this_transformation``
            this_transformation = sector.other.pop() # remove transformation from ``sector.other``

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
            for poly, hide_container in zip(sector.other, other_polynomials_regulator_hide_containers):
                hide_container.coeffs = poly.coeffs # coeffs had been hidden in `other_polynomials_name_hide_containers` --> do not overwrite
                decomposition.unhide(poly, hide_container)


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

            # apply ``this_transformation`` to ``this_sector_remainder_expression`` BEFORE taking derivatives
            for variable_index, integration_variable in enumerate(integration_variables):
                replacement = sp.sympify( Polynomial(this_transformation.expolist[variable_index:variable_index+1,:], [1], integration_variables) )
                this_sector_remainder_expression = this_primary_sector_remainder_expression.replace(variable_index, replacement)
            this_sector_remainder_expression = Expression(sp.sympify(this_sector_remainder_expression), symbols_remainder_expression)

            # define ``cal_I``, the part of the integrand that does not lead to poles
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

            # expand poles
            integrand_summands = []
            for i,(regular,singular) in enumerate(zip(regular_parts, pole_parts)):
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

                integrand_summands.append( Product(singular_expanded,regular_expanded,copy=False) )

            integrand = Sum(*integrand_summands, copy=False)

            # define the CFunctions for FORM
            # TODO: How to determine which derivatives of the user input ``functions`` are needed? How to communicate it to the user?
            all_functions = list(functions)
            functions_for_insertion = []

            # compute the required derivatives
            derivatives = {}
            def update_derivatives(basename, derivative_tracker, full_expression):
                derivatives[basename] = full_expression # include undifferentiated expression
                all_functions.append(basename) # define the symbol as CFunction in FORM
                functions_for_insertion.append(basename) # "mark" for insertion
                for multiindex,expression in derivative_tracker.compute_derivatives(full_expression).items():
                    name = _derivative_muliindex_to_name(basename, multiindex)
                    derivatives[name] = expression
                    all_functions.append(name) # define the symbol as CFunction in FORM
                    functions_for_insertion.append(name) # "mark" for insertion

            #  - for the contour deformation
            if contour_deformation_polynomial is not None:
                update_derivatives(
                    contour_deformation_polynomial, # basename
                    symbolic_contour_deformation_polynomial, # derivative tracker
                    sector.cast[contour_deformation_polynomial_index] # full expression
                )

            update_derivatives(basename=FORM_names['cal_I'], derivative_tracker=symbolic_cal_I, full_expression=cal_I)


            # generate the function definitions for the insertion in FORM
            FORM_function_definitions = ''.join(
                _make_FORM_function_definition(name, expression, limit=10**6)
                for name, expression in derivatives.items()
            )
            form_insertions = _make_FORM_list(derivatives.keys())

            # generate list over all occuring orders in the regulators
            regulator_powers = list( rangecomb(np.zeros_like(requested_orders), requested_orders+highest_poles_current_sector) )
            number_of_orders = len(regulator_powers)

            # generate the definitions of the FORM preprocessor variables "shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex'"
            regulator_powers = _make_FORM_shifted_orders(regulator_powers)

            # parse template file "sector.h"
            template_replacements['functions'] = _make_FORM_list(all_functions)
            template_replacements['function_definitions'] = FORM_function_definitions
            template_replacements['integrand'] = integrand
            template_replacements['integrand_container_initializer'] = _make_FORM_Series_initilization(-highest_poles_current_sector, requested_orders, sector_index)
            template_replacements['highest_regulator_poles'] = _make_FORM_list(highest_poles_current_sector)
            template_replacements['regulator_powers'] = regulator_powers
            template_replacements['number_of_orders'] = number_of_orders
            template_replacements['functions_for_insertion'] = _make_FORM_list(functions_for_insertion)
            parse_template_file(os.path.join(template_sources, 'name', 'codegen', 'sector.h'), # source
                                os.path.join(target_directory, name, 'codegen', 'sector%i.h' % sector_index), # dest
                                template_replacements)
            # TODO: adapt and parse "contour_deformation.h"

    # parse the template file "integrands.hpp"
    template_replacements['number_of_sectors'] = sector_index
    template_replacements['sector_includes'] = ''.join( '#include <%s/integrands/sector_%i.hpp>\n' % (name, i) for i in range(1,sector_index+1) )
    template_replacements['sectors_initializer'] = ','.join( 'integrand_of_sector_%i' % i for i in range(1,sector_index+1) )
    parse_template_file(os.path.join(template_sources, 'name', 'integrands', 'integrands.hpp'), # source
                        os.path.join(target_directory,  name , 'integrands', 'integrands.hpp'), # dest
                        template_replacements)
