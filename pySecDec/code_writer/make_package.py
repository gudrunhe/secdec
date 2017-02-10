"""
This module implements the main program - the function
:func:`.make_package`.

"""

from __future__ import print_function
from ..metadata import version, git_id
from ..misc import sympify_symbols, rangecomb
from ..algebra import _Expression, Expression, Polynomial, \
                      ExponentiatedPolynomial, Pow, Product, \
                      ProductRule, Function, Sum
from .. import decomposition
from ..matrix_sort import iterative_sort, Pak_sort
from ..subtraction import integrate_pole_part, integrate_by_parts, pole_structure as compute_pole_structure
from ..expansion import expand_singular, expand_Taylor, expand_sympy
from ..misc import lowest_order, det
from .template_parser import parse_template_file, parse_template_tree
from itertools import chain
from time import strftime
from re import match
import numpy as np
import sympy as sp
import sys, os

# The only public object this module provides is the function `make_package`.
# The module is organized in multiple sections dedicated to specific tasks
# to be addressed while writing the c++ package.
# To find the main function `make_package`, it is easiest to full-text search
# for "def make_package(".

_sympy_one = sp.sympify(1)

# define the internal names to be used in FORM
internal_prefix = 'SecDecInternal'
FORM_names = dict(
    cal_I=internal_prefix+'CalI',
    real_part=internal_prefix+'RealPart',
    cast_polynomial=internal_prefix+'PolynomialToDecompose',
    other_polynomial=internal_prefix+'OtherPolynomial',
    remainder_expression=internal_prefix+'RemainderExpression',
    error_token=internal_prefix+'ErrorToken',
    contourdef_transform=internal_prefix+'ContourdefDeformation',
    contourdef_Jacobian=internal_prefix+'ContourdefJacobian',
    deformed_variable=internal_prefix+'Deformed',
    deformation_parameter_i=internal_prefix+'Lambda%i',
    additional_deformation_factor=internal_prefix+'AdditionalDeformationFactor'
)

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
        if target_type is ExponentiatedPolynomial:
            if type(expression) == ExponentiatedPolynomial:
                expression.exponent = sp.sympify(str(expression.exponent))
                expression.coeffs = np.array([ sp.sympify(str(coeff)) for coeff in expression.coeffs ])
            elif type(expression) == Polynomial:
                expression = ExponentiatedPolynomial(expression.expolist,
                                                     np.array([ sp.sympify(str(coeff)) for coeff in expression.coeffs ]),
                                                     _sympy_one, # exponent
                                                     polysymbols, copy=False)
            else:
                expression = sp.sympify(str(expression))
                if expression.is_Pow:
                    assert len(expression.args) == 2
                    expression_base = Polynomial.from_expression(expression.args[0], polysymbols)
                    expression_exponent = sp.sympify(str(expression.args[1]))
                    expression = ExponentiatedPolynomial(expression_base.expolist,
                                                         expression_base.coeffs,
                                                         expression_exponent, polysymbols, copy=False)
                else:
                    expression = Polynomial.from_expression(expression, polysymbols)
                    expression = ExponentiatedPolynomial(expression.expolist,
                                                         expression.coeffs,
                                                         _sympy_one, # exponent
                                                         polysymbols, copy=False)
        elif target_type is _Expression:
            if not isinstance(expression, _Expression):
                expression, functions = Expression(expression, polysymbols, follow_functions=True)
        else:
            raise RuntimeError('`_parse_expressions` is only implemented for `target_type`s in %s, not for %s. (raised while parsing `%s`)' \
                                % (set([_Expression, ExponentiatedPolynomial]), target_type, name_of_make_argument_being_parsed))
        parsed.append(expression if target_type is ExponentiatedPolynomial else (expression,functions))
    return parsed

def _validate(name):
    '''
    Check validity of `name` for usage in FORM,
    c++, and the file system.
    Restrictions are as follows:
     o FORM: no underscores
     o FORM, c++: the first character must be
       in [A-Z, a-z]; i.e. no leading numbers
     o all: no special characters; i.e. only
       characters from the set [A-Z, a-z, 0-9]
     o `name` must not begin with the internal
       prefix

    '''
    if match(r'^[A-Z,a-z]+[A-Z,a-z,0-9]*$', name) is None:
        raise NameError('"%s" cannot be used as symbol' % name)
    if name.startswith(internal_prefix):
        raise NameError('Symbol names must not start with "%s"' % internal_prefix)

def _convert_input(name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   real_parameters, complex_parameters, form_optimization_level,
                   form_work_space, form_insertion_depth, contour_deformation_polynomial,
                   positive_polynomials, decomposition_method):
    'Get the data types right.'

    # parse symbols
    integration_variables = sympify_symbols(list(integration_variables), 'All `integration_variables` must be symbols.')
    regulators = sympify_symbols(list(regulators), 'All `regulators` must be symbols.')
    polynomial_names = sympify_symbols(list(polynomial_names), 'All `polynomial_names` must be symbols.')
    positive_polynomials = sympify_symbols(list(positive_polynomials), 'All `positive_polynomials` must be symbols.')
    real_parameters= sympify_symbols(list(real_parameters), 'All `real_parameters` must be symbols.')
    complex_parameters = sympify_symbols(list(complex_parameters), 'All `complex_parameters` must be symbols.')
    functions = list(functions); sympify_symbols(functions, 'All `functions` must be symbols.'); functions = set(str(f) for f in functions)
    if contour_deformation_polynomial is not None:
        contour_deformation_polynomial = sympify_symbols([contour_deformation_polynomial], '`contour_deformation_polynomial` must be a symbol.')[0]

    # check validity of symbol names and `name`
    for symbol in chain(integration_variables, regulators, polynomial_names, real_parameters, complex_parameters, functions,
                        [contour_deformation_polynomial] if contour_deformation_polynomial is not None else []):
        _validate( str(symbol) )

    # check that all the `positive_polynomials` also appear in `polynomial_names`
    for polyname in positive_polynomials:
        assert polyname in polynomial_names, '"%s" found in `positive_polynomials` but not in `polynomial_names` (%s)' % (polyname,polynomial_names)

    # define the symbols of the different classes of `_Expression`s
    symbols_polynomials_to_decompose = integration_variables + regulators + polynomial_names
    symbols_remainder_expression = integration_variables + regulators + polynomial_names
    all_symbols = symbols_other_polynomials = integration_variables + regulators + polynomial_names

    # check format of `requested_orders` --> must be 1d and have the length as `regulators`
    requested_orders = np.array(requested_orders)
    assert len(requested_orders.shape) == 1, '`requested_orders` has the wrong shape. (is %s; should be (%i,))' % (requested_orders.shape,len(regulators))
    assert requested_orders.shape[0] == len(regulators), 'The length of `requested_orders` (%i) must match the length of `regulators` (%i)' % (len(requested_orders),len(regulators))

    # parse expressions
    polynomials_to_decompose = _parse_expressions(polynomials_to_decompose, symbols_polynomials_to_decompose, ExponentiatedPolynomial, 'polynomials_to_decompose')
    other_polynomials = _parse_expressions(other_polynomials, symbols_other_polynomials, ExponentiatedPolynomial, 'other_polynomials')
    remainder_expression, function_calls = _parse_expressions([remainder_expression], symbols_remainder_expression, _Expression, 'remainder_expression')[0]

    # the exponents be polynomials in the regulators
    for poly in polynomials_to_decompose + other_polynomials:
        try:
            Polynomial.from_expression(poly.exponent, regulators)
        except Exception as error:
            error.args = tuple(['The exponents of the `polynomials_to_decompose` and the `other_polynomials` must be polynomials in the regulators. Error while checking: "%s"' % poly] + [arg for arg in error.args])
            raise
        exponent = Polynomial.from_expression(poly.exponent, integration_variables)
        assert (exponent.expolist == 0).all(), 'The exponents of the `polynomials_to_decompose` and the `other_polynomials` must not depend on the `integration_variables`. Error while checking: "%s"' % poly

    # convert ``prefactor`` to sympy expression
    prefactor = sp.sympify(prefactor)

    # convert ``requested_orders`` to numpy array
    requested_orders = np.array(requested_orders)

    # convert ``form_insertion_depth`` to integer and check nonnegativity
    assert form_insertion_depth == int(form_insertion_depth), '`form_insertion_depth` must be an integer.'
    assert form_insertion_depth >= 0, '`form_insertion_depth` must not be negative.'

    # check that neither the `remainder_expression` nor the `polynomials_to_decompose` refer to any of the `polynomial_names`
    str_error = internal_prefix + 'Error'
    error = sp.symbols(str_error)
    for _ in polynomial_names:
        remainder_expression = remainder_expression.replace(-1, error)
    str_remainder_expression = str(remainder_expression)
    if str_error in str_remainder_expression:
        raise ValueError('The `polynomial_names` %s cannot be used in the `remainder_expression`.' % polynomial_names)
    for poly in polynomials_to_decompose:
        for _ in polynomial_names:
            poly = poly.replace(-1, error)
        str_poly = str(poly)
        if str_error in str_poly:
            raise ValueError('The `polynomial_names` %s cannot be used in the `polynomials_to_decompose`.' % polynomial_names)

    return (name, integration_variables, regulators,
            requested_orders, polynomials_to_decompose, polynomial_names,
            other_polynomials, prefactor, remainder_expression, functions, function_calls,
            real_parameters, complex_parameters, form_optimization_level,
            form_work_space, form_insertion_depth, contour_deformation_polynomial, positive_polynomials,
            decomposition_method, symbols_polynomials_to_decompose, symbols_other_polynomials,
            symbols_remainder_expression, all_symbols)


# ---------------------------------- decomposition ----------------------------------
def get_decomposition_routines(name, normaliz, workdir):
    '''
    Return a dictionary with the functions
    performing the primary and the secondary
    decomposition.
    Along with the name, the path to the executable
    of normaliz and a temporary directory are passed
    to this function.

    '''
    _decomposition_strategies = dict(
                                        iterative=           dict(
                                                                     primary=decomposition.iterative.primary_decomposition,
                                                                     secondary=decomposition.iterative.iterative_decomposition
                                                                 ),
                                        geometric=           dict(
                                                                     primary=lambda sector, indices: [decomposition.geometric.Cheng_Wu(sector, indices[-1])],
                                                                     secondary=lambda sector, indices: decomposition.geometric.geometric_decomposition(sector, indices, normaliz, workdir)
                                                                 ),
                                        iterative_no_primary=dict(
                                                                     primary=lambda sector, indices: [sector], # no primary decomposition
                                                                     secondary=decomposition.iterative.iterative_decomposition
                                                                 )
                                    )
    return _decomposition_strategies[name]


# -------------------------------- template parsing ---------------------------------
def _parse_global_templates(name, regulators, polynomial_names,
                            real_parameters, complex_parameters, form_optimization_level,
                            form_work_space, form_insertion_depth, requested_orders,
                            contour_deformation_polynomial, nested_series_type,
                            enforce_complex):
    '''
    Create the `target_directory` (given by `name`) and return the two
    optional arguments passed to :func:`parse_template_tree`.

    '''
    # initialize template replacements
    template_replacements = dict(
                                     name = name,
                                     number_of_real_parameters = len(real_parameters),
                                     real_parameters = _make_FORM_list(real_parameters),
                                     names_of_real_parameters = _make_cpp_list(real_parameters),
                                     number_of_complex_parameters = len(complex_parameters),
                                     complex_parameters = _make_FORM_list(complex_parameters),
                                     names_of_complex_parameters = _make_cpp_list(complex_parameters),
                                     have_complex_parameters = len(complex_parameters) > 0,
                                     number_of_regulators = len(regulators),
                                     regulators = _make_FORM_list(regulators),
                                     names_of_regulators = _make_cpp_list(regulators),
                                     polynomial_names = _make_FORM_list(regulators),
                                     form_optimization_level = form_optimization_level,
                                     form_work_space = form_work_space,
                                     form_insertion_depth = form_insertion_depth,
                                     contour_deformation = int(contour_deformation_polynomial is not None),
                                     requested_orders = _make_FORM_list(requested_orders),
                                     nested_series_type = nested_series_type,
                                     pySecDec_version = version,
                                     python_version = sys.version,
                                     pySecDec_git_id = git_id,
                                     date_time = strftime("%a %d %b %Y %H:%M"),
                                     enforce_complex_return_type=int(bool(enforce_complex)) # make sure that this is either ``0`` or ``1``
                                )

    # configure template parser
    file_renamings = {
                          # replace "name" by the name of the integral
                          'name' : name,
                          'name.cpp' : name + '.cpp',

                          # the files below are specific for each sector --> do not parse globally
                          'contour_deformation.h' : None,
                          'sector.h' : None,

                          # "integrands.cpp", "name.hpp", "prefactor.cpp", "pole_structures.cpp", and "functions.hpp" can only be written after the decomposition is completed
                          'integrands.cpp' : None,
                          'name.hpp' : None,
                          'prefactor.cpp' : None,
                          'pole_structures.cpp' : None,
                          'functions.hpp' : None
                     }

    # the files below are only relevant for contour deformation --> do not parse if deactivated
    if contour_deformation_polynomial is None:
        for filename in ['contour_deformation.h','write_contour_deformation.frm']:
            file_renamings[filename] = None

    # get path to the directory with the template files (path relative to directory with this file: "./templates/")
    template_sources = os.path.join(os.path.split(os.path.abspath(__file__))[0],'templates')

    # initialize the target directory with the sector independent files
    parse_template_tree(template_sources, name, template_replacements, file_renamings)

    # return parser options
    return template_sources, template_replacements, file_renamings


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

def _make_FORM_function_definition(name, expression, args, limit):
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

    :param args:
        iterable of sympy symbols or None;
        The arguments of the function.

    :param limit:
        integer;
        The maximum number of characters in each
        subexpression.

    '''
    if args is None:
        FORM_args_left_hand_side = FORM_args_right_hand_side = ''
    else:
        FORM_args_left_hand_side = '(' + _make_FORM_list(str(arg)+'?' for arg in args) + ')'
        FORM_args_right_hand_side = '(' + _make_FORM_list(str(arg) for arg in args) + ')'
    codelines = []

    def recursion(name, expression):
        str_expression = str(expression)

        if len(str_expression) <= limit: # no need to split
            codelines.append( "  Id %s = %s;" % (name+FORM_args_left_hand_side,str_expression) )
            return

        if type(expression) == ProductRule:
            expression = expression.to_sum()

        if type(expression) == Sum:
            name_next_level = internal_prefix+'fDUMMY'+name+'Part'
            codelines.append(
                "  Id %s = %s;" % (
                    name+FORM_args_left_hand_side,
                    '+'.join( name_next_level+str(i)+FORM_args_right_hand_side for i in range(len(expression.summands)) )
                )
            )
            for i,summand in enumerate(expression.summands):
                recursion(name_next_level+str(i), summand)
            return

        if type(expression) == Product:
            name_next_level = internal_prefix+'fDUMMY'+name+'Part'
            codelines.append(
                "  Id %s = %s;" % (
                    name+FORM_args_left_hand_side,
                    '*'.join( name_next_level+str(i)+FORM_args_right_hand_side for i in range(len(expression.factors)) )
                )
            )
            for i,factor in enumerate(expression.factors):
                recursion(name_next_level+str(i), factor)
            return

        # rescue: print warning and write unsplit expression
        print( 'WARNING: Could not split "%s" (not implemented for %s)' % (name,type(expression)) )
        codelines.append( "  Id %s = %s;" % (name+FORM_args_left_hand_side,str_expression) )
        return

    recursion(name, expression)
    codelines.append('') # empty line
    return (  "\n".join(codelines)  ).replace('**','^')

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
    :meth:`pySecDec.algebra.Function.compute_derivatives`
    to the form ``d...d<basename>d<index>...d<index>``.

    Example:

    >>> _derivative_muliindex_to_name('f', (1,2,1))
    ddddfd0d1d1d2

    '''
    prefix = 'd' * sum(multiindex)
    suffix = ''.join(('d' + str(i)) * depth for i,depth in enumerate(multiindex))
    return prefix + basename + suffix


# ---------------------------------- write c++ code ---------------------------------
def _make_CXX_Series_initialization(regulator_names, min_orders, max_orders, sector_ID, contour_deformation):
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
            outstr_tail = '},true,#@%sDblquote@#%s#@%sDblquote@#}' % (internal_prefix,regulator_names[regulator_index],internal_prefix)
            return ''.join( (outstr_head, ','.join(outstr_body_snippets), outstr_tail) )
        else: # regulator_index == last_regulator_index; i.e. processing last regulator
            outstr_head = '{%i,%i,{{' % (min_orders[regulator_index],max_orders[regulator_index])
            outstr_body_snippets = []
            for this_regulator_order in range(min_orders[regulator_index],max_orders[regulator_index]+1):
                current_orders[regulator_index] = this_regulator_order
                if contour_deformation:
                    outstr_body_snippets.append(
                        '''%(sector_ID)i,\{%(order)s\},sector_%(sector_ID)i_order_%(cpp_order)s_numIV,sector_%(sector_ID)i_order_%(cpp_order)s_integrand,
                           sector_%(sector_ID)i_order_%(cpp_order)s_contour_deformation_polynomial,
                           sector_%(sector_ID)i_order_%(cpp_order)s_maximal_allowed_deformation_parameters''' \
                        % dict(sector_ID=sector_ID,cpp_order=multiindex_to_cpp_order(current_orders),order=_make_FORM_list(current_orders))
                    )
                else:
                    outstr_body_snippets.append(
                        '%(sector_ID)i,\{%(order)s\},sector_%(sector_ID)i_order_%(cpp_order)s_numIV,sector_%(sector_ID)i_order_%(cpp_order)s_integrand' \
                        % dict(sector_ID=sector_ID,cpp_order=multiindex_to_cpp_order(current_orders),order=_make_FORM_list(current_orders))
                    )
            outstr_tail = '}},true,#@%sDblquote@#%s#@%sDblquote@#}' % (internal_prefix,regulator_names[regulator_index],internal_prefix)
            return ''.join( (outstr_head, '},{'.join(outstr_body_snippets), outstr_tail) )

    return recursion(0)

def _make_cpp_list(python_list):
    '''
    Convert a python list to a string to be used
    in c++ initializer list for a vector.

    Example: ``['a', 'b', 'c'] --> '"a","b","c"'``

    '''
    joined_inner_part = '","'.join(str(item) for item in python_list)
    if len(joined_inner_part) == 0:
        return ''
    else:
        return '"' + joined_inner_part + '"'

def _make_prefactor_function(expanded_prefactor, real_parameters, complex_parameters):
    regulators = expanded_prefactor.polysymbols
    last_regulator_index = len(regulators) - 1

    def recursion(regulator_index,expression):
        orders = expression.expolist[:,regulator_index]
        min_order = orders.min()
        max_order = orders.max()
        if regulator_index < last_regulator_index:
            outstr_body_snippets = []
            outstr_head = '{%i,%i,{' % (min_order,max_order)
            for coeff in expression.coeffs:
                outstr_body_snippets.append( recursion(regulator_index + 1, coeff) )
            outstr_tail = '},true' if expression.truncated else '},false'
            outstr_tail += ',"%s"}' % regulators[regulator_index]
            return ''.join( (outstr_head, ','.join(outstr_body_snippets), outstr_tail) )
        else: # regulator_index == last_regulator_index; i.e. processing last regulator
            outstr_head = '{%i,%i,{{' % (min_order,max_order)
            outstr_body_snippets = []
            for coeff in expression.coeffs:
                outstr_body_snippets.append( str(coeff.evalf(20)) )
            outstr_tail = '}},true' if expression.truncated else '}},false'
            outstr_tail += ',"%s"}' % regulators[regulator_index]
            return ''.join( (outstr_head, '},{'.join(outstr_body_snippets), outstr_tail) )

    parameter_definitions = []
    parameter_undefs = []
    for real_or_complex in ('real','complex'):
        for i,p in enumerate(eval(real_or_complex+'_parameters')):
            parameter_definitions.append('#define %s %s_parameters.at(%i)' % (p,real_or_complex,i))
            parameter_undefs.append('#undef %s' % p)
    parameter_definitions = '\n'.join(parameter_definitions)
    parameter_undefs = '\n'.join(parameter_undefs)

    code = parameter_definitions + '\nreturn ' + recursion(0,expanded_prefactor) + ';\n' + parameter_undefs
    return code.replace('\n', '\n        ')

def _make_CXX_function_declaration(function_name, number_of_arguments):
    '''
    Write the declaration of a c++ function with
    name `function_name` and `number_of_arguments`
    arguments. The function's return type is set
    to "integrand_return_t". The argument types
    are templated.

    '''
    if number_of_arguments == 0:
        return '    integrand_return_t ' + function_name + '();\n'

    template_arguments = ', '.join('typename T%i' % i for i in range(number_of_arguments))
    arguments = ', '.join('T%i arg%i' % (i,i) for i in range(number_of_arguments))
    return '    template<' + template_arguments + '>\n    integrand_return_t ' + function_name + '(' + arguments + ');\n'


# --------------------------------- algebra helper ---------------------------------
class RealPartFunction(Function):
    '''
    Symbolic function that takes exactly one argument
    and with the additional property that the
    derivative commutes to the inside.

    '''
    def __init__(self, symbol, *arguments, **kwargs):
        assert len(arguments) == 1, 'A `RealPartFunction` takes exactly one argument'
        super(RealPartFunction, self).__init__(symbol, *arguments, **kwargs)

    def derive(self, index):
        return RealPartFunction(self.symbol, self.arguments[0].derive(index), copy=False)


# ---------------------------------- main function ----------------------------------
def make_package(name, integration_variables, regulators, requested_orders,
                 polynomials_to_decompose, polynomial_names=[], other_polynomials=[],
                 prefactor=1, remainder_expression=1, functions=[], real_parameters=[],
                 complex_parameters=[], form_optimization_level=2, form_work_space='500M',
                 form_insertion_depth=5, contour_deformation_polynomial=None, positive_polynomials=[],
                 decomposition_method='iterative_no_primary', normaliz_executable='normaliz',
                 normaliz_workdir='normaliz_tmp', enforce_complex=False, split=False, ibp_power_goal=-1):
    r'''
    Decompose, subtract and expand an expression.
    Return it as c++ package.

    .. seealso::
        In order to decompose a loop integral,
        use the function
        :func:`pySecDec.loop_integral.loop_package`.

    :param name:
        string;
        The name of the c++ namepace and the output
        directory.

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

    :param functions:
        iterable of strings or sympy symbols, optional;
        Function symbols occuring in `remainder_expression`,
        e.g.``['f']``.

        .. note::
            The power function `pow` and the logarithm
            `log` are already defined by default. The
            `log` uses the nonstandard continuation
            from a negative imaginary part on the negative
            real axis (e.g. ``log(-1) = -i*pi``).

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
        How deep FORM should try to resolve nested function
        calls. Default: ``5``.

    :param contour_deformation_polynomial:
        string or sympy symbol, optional;
        The name of the polynomial in `polynomial_names`
        that is to be continued to the complex plane
        according to a :math:`- i\delta` prescription.
        For loop integrals, this is the second Symanzik
        polynomial ``F``.
        If not provided, no code for contour deformation
        is created.

    :param positive_polynomials:
        iterable of strings or sympy symbols, optional;
        The names of the polynomials in `polynomial_names`
        that should always have a positive real part.
        For loop integrals, this applies to the first Symanzik
        polynomial ``U``.
        If not provided, no polynomial is checked for
        positiveness.
        If `contour_deformation_polynomial` is ``None``, this
        parameter is ignored.

    :param decomposition_method:
        string, optional;
        The strategy to decompose the polynomials. The
        following strategies are available:

        * 'iterative_no_primary' (default)
        * 'iterative'
        * 'geometric'

        'iterative' and 'geometric' are only valid for
        loop integrals. An end user should always use
        the default 'iterative_no_primary' here.
        In order to compute loop integrals, please use the
        function :func:`pySecDec.loop_integral.loop_package`.

    :param normaliz_executable:
        string, optional;
        The command to run `normaliz`. `normaliz` is only
        required if `decomposition_method` is set to
        'geometric'.
        Default: 'normaliz'

    :param normaliz_workdir:
        string, optional;
        The working directory for `normaliz`. This directory
        is automatically created and deleted again after
        `normaliz` finishes.
        Default: 'normaliz_tmp'

    :param enforce_complex:
        bool, optional;
        Whether or not the generated integrand functions
        should have a complex return type even though
        they might be purely real.
        The return type of the integrands is automatically
        complex if `contour_deformation` is ``True`` or
        if there are `complex_parameters`. In other cases,
        the calculation can typically be kept purely real.
        Most commonly, this flag is needed if
        ``log(<negative real>)`` occurs in one of the
        integrand functions. However, `pySecDec` will suggest
        setting this flag to ``True`` in that case.
        Default: ``False``

    :param split:
        bool or integer, optional;
        Whether or not to split the integration domain
        in order to map singularities from :math:`1` to
        :math:`0`. Set this option to ``True`` if you have
        singularties when one or more integration variables
        are one. If an integer is passed, that integer is
        used as seed to generate the splitting point.
        Default: ``False``

    :param ibp_power_goal:
        integer, optional;
        The `power_goal` that is forwarded to
        :func:`.integrate_by_parts`.

        This option controls how the subtraction terms are
        generated. Setting it to ``-numpy.inf`` disables
        :func:`.integrate_by_parts`, while ``0`` disables
        :func:`.integrate_pole_part`.

        .. seealso::
            To generate the subtraction terms, this function
            first calls :func:`.integrate_by_parts` for each
            integration variable with the give `ibp_power_goal`.
            Then :func:`.integrate_pole_part` is called.

        Default: ``-1``

    '''
    print('running "make_package" for "' + name + '"')

    # convert input data types to the data types we need
    name, integration_variables, regulators, \
    requested_orders, polynomials_to_decompose, polynomial_names, \
    other_polynomials, prefactor, remainder_expression, functions, \
    function_calls, real_parameters, complex_parameters, \
    form_optimization_level, form_work_space, form_insertion_depth, \
    contour_deformation_polynomial, positive_polynomials, decomposition_method, \
    symbols_polynomials_to_decompose, symbols_other_polynomials, \
    symbols_remainder_expression, all_symbols = \
    _convert_input(name, integration_variables, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   real_parameters, complex_parameters, form_optimization_level,
                   form_work_space, form_insertion_depth, contour_deformation_polynomial,
                   positive_polynomials, decomposition_method)

    # construct the c++ type "nested_series_t"
    # for two regulators, the resulting code should read:
    # "secdecutil::Series<secdecutil::Series<T>>"
    nested_series_type = 'secdecutil::Series<' * len(regulators) + 'T' + '>' * len(regulators)

    # configure the template parser and parse global files
    template_sources, template_replacements, file_renamings = \
        _parse_global_templates(
        name, regulators, polynomial_names,
        real_parameters, complex_parameters, form_optimization_level,
        form_work_space, form_insertion_depth, requested_orders,
        contour_deformation_polynomial, nested_series_type,
        enforce_complex
    )

    # get the highest poles from the ``prefactor``
    highest_prefactor_pole_orders = -np.array([lowest_order(prefactor, regulator) for regulator in regulators])

    # compute the required expansion order accounting for the prefactor
    required_orders = requested_orders + highest_prefactor_pole_orders

    # get the decomposition routines
    strategy = get_decomposition_routines(decomposition_method, normaliz_executable, normaliz_workdir)

    # define the monomials "x0", "x1", "x2", ... to keep track of the transformations
    one = Polynomial([[0]*len(symbols_other_polynomials)], [1], symbols_other_polynomials)
    transformations = []
    for i in range(len(integration_variables)):
        to_append = one.copy()
        to_append.expolist[:,i] = 1
        transformations.append(to_append)

    # define an error token that is multiplied to expressions that should evaluate to zero
    error_token = sp.symbols(FORM_names['error_token'])

    # make a copy of the `integration_variables` for later reference
    all_integration_variables = list(integration_variables)

    # intialize the c++ declarations of the `functions`
    function_declarations = set()

    have_dummy_functions = True if functions else False

    # initialize the decomposition
    initial_sector = decomposition.Sector(polynomials_to_decompose, other_polynomials + transformations)

    # if splitting desired, implement it as additional primary decomposition
    if split:
        # cannot split when using the geometric decomposition method because the integration interval is [0,inf] after the primary decomposition
        if decomposition_method == 'geometric':
            raise ValueError('Cannot have ``split=True`` and ``decomposition_method="geometric"``. You probably want to try ``split=True`` and ``decomposition_method="iterative"``')

        original_decomposition_strategies = strategy

        def primary_decomposition_with_splitting(sector, indices):
            # investigate symmetries before the split
            if use_symmetries:
                primary_sectors = list(  original_decomposition_strategies['primary'](sector, indices)  )
                print('number of primary sectors before investigating symmetries:', len(primary_sectors))
                primary_sectors = decomposition.squash_symmetry_redundant_sectors(primary_sectors, iterative_sort)
                primary_sectors = decomposition.squash_symmetry_redundant_sectors(primary_sectors, Pak_sort)
                print('number of primary sectors after investigating symmetries:', len(primary_sectors))
            else:
                primary_sectors = original_decomposition_strategies['primary'](sector, indices)
            for output_sector in primary_sectors:
                yield output_sector

        def secondary_decomposition_with_splitting(sector, indices):
            # split and decompose the `sector`
            for split_sector in decomposition.splitting.split_singular(sector, split, indices):
                for decomposed_sector in original_decomposition_strategies['secondary'](split_sector, indices):
                    yield decomposed_sector

        # apply modifications to the decomposition strategy
        strategy = dict(primary=primary_decomposition_with_splitting, secondary=secondary_decomposition_with_splitting)

    # initialize the counter
    sector_index = 0

    # initialize the global lowest orders
    lowest_orders = requested_orders.copy()

    # initialize list of pole structures
    pole_structures = []

    # define the imaginary unit
    imaginary_unit = sp.sympify('I')

    # define the dummy names for the `other_polynomials`
    names_other_polynomials = [FORM_names['other_polynomial'] + str(i) for i in range(len(other_polynomials))]

    # we must backwards traverse the `polynomial_names` --> create the reversed list once and for all
    reversed_polynomial_names = list(polynomial_names) # copy
    reversed_polynomial_names.reverse()

    if contour_deformation_polynomial is not None:
        # get the index of the `contour_deformation_polynomial` in `polynomial_names`
        str_contour_deformation_polynomial = str(contour_deformation_polynomial)
        contour_deformation_polynomial_index = 0
        while str(polynomial_names[contour_deformation_polynomial_index]) != str_contour_deformation_polynomial:
            contour_deformation_polynomial_index += 1
            if len(polynomial_names) <= contour_deformation_polynomial_index:
                raise IndexError('Could not find the `contour_deformation_polynomial` "%s" in `polynomial_names`.' % str_contour_deformation_polynomial)

    def parse_exponents_and_coeffs(sector, symbols_polynomials_to_decompose, symbols_other_polynomials):
        #  - in ``sector.cast``
        for product in sector.cast:
            mono, poly = product.factors
            mono.exponent = Polynomial.from_expression(mono.exponent, symbols_polynomials_to_decompose)
            poly.exponent = Polynomial.from_expression(poly.exponent, symbols_polynomials_to_decompose)
            mono.coeffs = np.array([Expression(coeff, symbols_polynomials_to_decompose) for coeff in mono.coeffs])
            poly.coeffs = np.array([Expression(coeff, symbols_polynomials_to_decompose) for coeff in poly.coeffs])

        #  - in ``sector.other``
        for poly in sector.other:
            try:
                poly.exponent = Polynomial.from_expression(poly.exponent, symbols_other_polynomials)
            except AttributeError:
                pass # not an exponentiated polynomial --> nothing to do
            poly.coeffs = np.array([Expression(coeff, symbols_polynomials_to_decompose) for coeff in poly.coeffs])

    # investigate if we can take advantage of sector symmetries
    # We can simplify due to symmetries if the `remainder_expression`
    # does explicitly not depend on any integration variable and if
    # we do not split the integration region.
    remainder_expression_is_trivial = True
    for i,_ in enumerate(integration_variables):
        derivative = remainder_expression.derive(i).simplify()
        if not ( type(derivative) is Polynomial and np.array_equal(derivative.coeffs, [0]) and (derivative.expolist == 0).all() ):
            remainder_expression_is_trivial = False
            break
    use_symmetries = remainder_expression_is_trivial

    # Check that either the primary decomposition or the `remainder_expression` is trivial.
    # Note that the primary decomposition is specialized for loop integrals.
    if not remainder_expression_is_trivial and decomposition_method != 'iterative_no_primary':
        raise NotImplementedError('The primary decomposition is only implemented for loop integrals. Please perform the primary decomposition yourself and choose ``decomposition_method="iterative_no_primary"``.')

    if use_symmetries:
        # can investigate sector symmetries
        # we do not need `transformations` in that case --> remove from `initial_sector`
        for i in range(len(integration_variables)):
            initial_sector.other.pop()

    # symmetries are applied elsewhere if we split
    if use_symmetries and not split:
        # run primary decomposition and squash symmetry-equal sectors (using both implemented strategies)
        primary_sectors = list( strategy['primary'](initial_sector, range(len(integration_variables))) )
        print('number of primary sectors before investigating symmetries:', len(primary_sectors))
        primary_sectors = decomposition.squash_symmetry_redundant_sectors(primary_sectors, iterative_sort)
        primary_sectors = decomposition.squash_symmetry_redundant_sectors(primary_sectors, Pak_sort)
        print('number of primary sectors after investigating symmetries:', len(primary_sectors))

        # rename the `integration_variables` in all `primary_sectors` --> must have the same names in all primary sectors
        symbols_primary_sectors = primary_sectors[0].Jacobian.polysymbols
        for sector in primary_sectors:
            sector.Jacobian.polysymbols = list(symbols_primary_sectors) # copy
            for prod in sector.cast:
                for factor in prod.factors:
                    factor.polysymbols = list(symbols_primary_sectors) # copy
            for poly in sector.other:
                poly.polysymbols = list(symbols_primary_sectors) # copy

        # give one primary sector as representative for the global initialization
        primary_sectors_to_consider = [primary_sectors[0]]

    else: # if we cannot take advantage of symmetries
        primary_sectors_to_consider = strategy['primary'](initial_sector, range(len(integration_variables)))

    for primary_sector_index, primary_sector in enumerate(primary_sectors_to_consider):

        # primary decomposition removes one integration parameter --> redefine `integration_variables` and the symbols of the different classes of `_Expression`s
        integration_variables = primary_sector.Jacobian.polysymbols[:-len(regulators)-len(polynomial_names)]
        symbols_polynomials_to_decompose = symbols_other_polynomials = integration_variables + regulators
        symbols_remainder_expression = integration_variables + regulators
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

        # we later need ``0`` and ``1`` packed into specific types
        polynomial_zero = Polynomial(np.zeros([1,len(symbols_other_polynomials)], dtype=int), np.array([0]), symbols_other_polynomials, copy=False)
        polynomial_one = Polynomial(np.zeros([1,len(symbols_other_polynomials)], dtype=int), np.array([1]), symbols_other_polynomials, copy=False)
        pole_part_initializer = Pow(polynomial_one, -polynomial_one)

        if contour_deformation_polynomial is not None:
            symbolic_deformed_variable_names = [FORM_names['deformed_variable'] + str(original_varname) for original_varname in integration_variables]
            symbolic_deformation_factor_names = [FORM_names['additional_deformation_factor'] + str(original_varname) for original_varname in integration_variables]

            # Need all first and second derivatives of the `contour_deformation_polynomial`.
            # Since the `contour_deformation_polynomial` is left symbolic they are equal for every subsector after primary decomposition.
            symbolic_contour_deformation_polynomial = Function(str(contour_deformation_polynomial), *elementary_monomials)

            # compute the deformation of the integration parameters and its Jacobian matrix (see e.g. section 3.2 in arXiv:1601.03982):
            # ``z_k({x_k}) = x_k * (1 - i * lambda_k * (1-x_k) * Re(dF_dx_k))``, where "dF_dx_k" denotes the derivative of ``F`` by ``x_k``
            # Remark: The determinant of the Jacobian matrix is calculated numerically.
            deformation_parameters = [sp.symbols(FORM_names['deformation_parameter_i'] % i) for i in range(len(integration_variables))]
            deformed_integration_parameters = [
                                                     Sum(
                                                         elementary_monomials[k],
                                                         Product(
                                                            -imaginary_unit * deformation_parameters[k] * elementary_monomials[k],
                                                            1 - elementary_monomials[k],
                                                            RealPartFunction(FORM_names['real_part'], symbolic_contour_deformation_polynomial.derive(k), copy=False),
                                                         copy=False),
                                                     copy=False)
                                                     for k in range(len(integration_variables))
                                              ]

            # define symbols for ``z_k({x_k}) / x_k``
            deformation_factors = [
                                         Sum(
                                             polynomial_one,
                                             Product(
                                                -imaginary_unit * deformation_parameters[k] * polynomial_one,
                                                1 - elementary_monomials[k],
                                                RealPartFunction(FORM_names['real_part'], symbolic_contour_deformation_polynomial.derive(k), copy=False),
                                             copy=False),
                                         copy=False)
                                         for k in range(len(integration_variables))
                                  ]

            # define the `symbolic_deformed_variables` and the `symbolic_deformation_factors` to be inserted FORM
            symbolic_deformed_variables = [ Function(deformed_name, *elementary_monomials[:len(integration_variables)]) for deformed_name in symbolic_deformed_variable_names ]
            symbolic_deformed_variables.extend( (regulator for regulator in elementary_monomials[len(integration_variables):]) )
            symbolic_deformation_factors = [ Function(deformed_name, *elementary_monomials[:len(integration_variables)]) for deformed_name in symbolic_deformation_factor_names ]

            # generate the Jacobian determinant
            print('computing Jacobian determinant for primary sector', primary_sector_index)
            contourdef_Jacobian = np.empty((len(integration_variables), len(integration_variables)), dtype=object)
            for i in range(len(integration_variables)):
                for j in range(len(integration_variables)):
                    contourdef_Jacobian[i,j] = symbolic_deformed_variables[i].simplify().derive(j)
            contourdef_Jacobian_determinant = det(contourdef_Jacobian)

        # remove `polynomial_names` from the `remainder_expression`
        this_primary_sector_remainder_expression = remainder_expression
        for poly_name in reversed_polynomial_names:
            this_primary_sector_remainder_expression = this_primary_sector_remainder_expression.replace(-1, poly_name(*symbols_polynomials_to_decompose), remove=True)

        # If there is a nontrivial primary decomposition, remove the integration variables from the `remainder_expression`
        for i,var in enumerate(all_integration_variables):
            if var not in integration_variables:
                this_primary_sector_remainder_expression = this_primary_sector_remainder_expression.replace(i,1,remove=True)
                break

        # define symbols for the `polynomials_to_decompose` and the `other_polynomials --> shorter expressions and faster in python
        symbolic_other_polynomials = [
                                         Pow(
                                             Function(
                                                         FORM_names['other_polynomial'] + str(i),
                                                         *(elementary_monomials if contour_deformation_polynomial is None else symbolic_deformed_variables)
                                             ),
                                             Expression(primary_sector.other[i].exponent, symbols_polynomials_to_decompose),
                                             copy = False
                                         )
                                         for i in range(len(other_polynomials))
                                     ]
        symbolic_polynomials_to_decompose = []
        names_polynomials_to_decompose = []
        for i in range(len(polynomials_to_decompose)):
            try:
                poly_name = str(polynomial_names[i])
            except IndexError:
                poly_name = FORM_names['cast_polynomial'] + str(i)
            symbolic_polynomials_to_decompose.append(
                Pow(
                    Function(poly_name, *(elementary_monomials if contour_deformation_polynomial is None else symbolic_deformed_variables)),
                    Expression(primary_sector.cast[i].factors[1].exponent, symbols_polynomials_to_decompose),
                    copy = False
                )
            )
            names_polynomials_to_decompose.append(poly_name)

        if use_symmetries and not split:
            # search for symmetries throughout the secondary decomposition
            secondary_sectors = []
            for primary_sector in primary_sectors:
                secondary_sectors.extend( strategy['secondary'](primary_sector, range(len(integration_variables))) )
            print('total number of sectors before investigating symmetries:', len(secondary_sectors))
            # find symmetries using both implemented strategies
            secondary_sectors = decomposition.squash_symmetry_redundant_sectors(secondary_sectors, iterative_sort)
            secondary_sectors = decomposition.squash_symmetry_redundant_sectors(secondary_sectors, Pak_sort)
            print('total number of sectors after investigating symmetries', len(secondary_sectors))
        else:
            secondary_sectors = strategy['secondary'](primary_sector, range(len(integration_variables)))

        for sector in secondary_sectors:
            sector_index += 1

            print('writing FORM files for sector', sector_index)

            if not use_symmetries:
                # If we do not use symmetries, we have to take care of `transformations`
                # extract ``this_transformations``
                this_transformations = sector.other[-len(all_integration_variables):]
                # remove transformation from ``sector.other``
                sector.other = sector.other[:-len(all_integration_variables)]

            # insert the `polynomials_to_decompose`:
            # explicitly insert monomial part as ``<symbol> --> xi**pow_i * <symbol>``
            # ``[:,:-len(regulators)-len(polynomial_names)]`` is the part of the expolist that contains only the powers of the `integration_variables`;
            #     i.e. it excludes the powers of the `regulators` and the `polynomial_names` from the modification
            for poly in sector.other:
                for i,(polyname,prod) in enumerate(zip(polynomial_names,sector.cast)):
                    mono,_ = prod.factors
                    poly.expolist[:,:-len(regulators)-len(polynomial_names)] += \
                        np.einsum('i,k->ik', poly.expolist[:,-len(polynomial_names):][:,i], mono.expolist[0,:-len(regulators)-len(polynomial_names)])

            # remove `polynomial_names` - keep polynomial part symbolic as dummy function:
            #  - from `other_polynomials`
            for i in range(len(other_polynomials)):
                poly = sector.other[i]
                for poly_name in reversed_polynomial_names:
                    poly = poly.replace(-1, poly_name(*symbols_other_polynomials), remove=True)
                sector.other[i] = poly

            #  - from `polynomials_to_decompose`
            for i in range(len(polynomials_to_decompose)):
                # no dependence here, just remove the symbols
                prod = sector.cast[i]
                for poly_name in reversed_polynomial_names:
                    prod = prod.replace(-1, poly_name(*symbols_polynomials_to_decompose), remove=True)
                sector.cast[i] = prod

            #  - from `Jacobian`
            Jacobian = sector.Jacobian
            for poly_name in reversed_polynomial_names:
                Jacobian = Jacobian.replace(-1, poly_name(*symbols_other_polynomials), remove=True)

            #  - from `this_transformations`
            if not use_symmetries:
                for i in range(len(all_integration_variables)):
                    for poly_name in reversed_polynomial_names:
                        this_transformations[i] = this_transformations[i].replace(-1, poly_name(*symbols_other_polynomials), remove=True)

            # convert all exponents and coefficients to pySecDec expressions
            parse_exponents_and_coeffs(sector, symbols_polynomials_to_decompose, symbols_other_polynomials)

            # factorize
            #  - the polynomials in ``sector.other``
            for i,to_factorize in enumerate(sector.other):
                to_factorize = sector.other[i] = \
                    Product(
                        ExponentiatedPolynomial(
                            np.zeros([1,len(symbols_other_polynomials)], dtype=int), np.array([1]), to_factorize.exponent.copy(), symbols_other_polynomials, copy=False
                        ),
                        to_factorize
                    )
                decomposition.refactorize(to_factorize)

            #  - the Jacobian
            Jacobian = \
                    Product(
                        ExponentiatedPolynomial(
                            np.zeros([1,len(Jacobian.polysymbols)], dtype=int), np.array([1]), exponent=polynomial_one, polysymbols=Jacobian.polysymbols, copy=False
                        ),
                        ExponentiatedPolynomial(Jacobian.expolist, Jacobian.coeffs, polysymbols=Jacobian.polysymbols, exponent=polynomial_one, copy=False)
                    )
            decomposition.refactorize(Jacobian)

            # Apply ``this_transformation`` and the contour deformaion (if applicable) to
            # the `remainder_expression` BEFORE taking derivatives.
            # Introduce a symbol for the `remainder_expression` and insert in FORM.
            # Note: ``elementary_monomials[len(integration_variables):]`` are the regulators
            if remainder_expression_is_trivial:
                # `remainder_expression` does not depend on the integration variables in that case
                symbolic_remainder_expression_arguments = [polynomial_zero] * len(integration_variables) + elementary_monomials[len(integration_variables):]
            else:
                symbolic_remainder_expression_arguments = this_transformations + elementary_monomials[len(integration_variables):]
            symbolic_remainder_expression = Function(FORM_names['remainder_expression'], *symbolic_remainder_expression_arguments)

            # initialize the product of monomials for the subtraction
            monomial_factors = list(chain([Jacobian.factors[0]], (prod.factors[0] for prod in sector.cast), (prod.factors[0] for prod in sector.other)))
            monomials = Product(*monomial_factors, copy=False)

            # compute the pole structure
            pole_structures.append(  compute_pole_structure(monomials, *integration_variable_indices)  )

            if contour_deformation_polynomial is not None:
                # Apply the deformation ``z_k({x_k}) = x_k - i * lambda_k * x_k * (1-x_k) * Re(dF_dx_k)`` to the monomials.
                # Split as ``z_k({x_k}) = monomials[k] * <something in remainder_expression>``
                # where ``monomials[k] = x_k``; i.e. unchanged
                # and where ``<something in remainder_expression> z_k({x_k}) / x_k = 1 - i * lambda_k * (1-x_k) * Re(dF_dx_k)``.
                # That amounts to multiplying ``<something in remainder_expression> ** <exponent_of_monomial>`` to `remainder_expression`
                additional_deformation_factors = []
                monomial_powers = np.zeros(len(integration_variables), dtype=object)
                for factor in monomial_factors:
                    for k,exponent in enumerate(factor.expolist[0,:len(integration_variables)]):
                        monomial_powers[k] += factor.exponent*exponent
                for k,(deformation_factor,power) in enumerate(zip(symbolic_deformation_factors,monomial_powers)):
                    additional_deformation_factors.append \
                    (
                        Pow(
                               deformation_factor,
                               power
                           )
                    )
                additional_deformation_factor = Product(*additional_deformation_factors, copy=False)

            # define ``cal_I``, the part of the integrand that does not lead to poles
            # use the derivative tracking dummy functions for the polynomials --> faster
            cal_I = Product(symbolic_remainder_expression, *chain([Jacobian.factors[1]], symbolic_polynomials_to_decompose, symbolic_other_polynomials), copy=False)

            # multiply Jacobian determinant to `cal_I`
            if contour_deformation_polynomial is not None:
                symbolic_contourdef_Jacobian = Function(FORM_names['contourdef_Jacobian'], *elementary_monomials[:len(integration_variables)])
                symbolic_additional_deformation_factor = Function(FORM_names['additional_deformation_factor'], *elementary_monomials)
                cal_I = Product(symbolic_contourdef_Jacobian, symbolic_additional_deformation_factor, cal_I)

            # it is faster to use a dummy function for ``cal_I`` and substitute back in FORM
            symbolic_cal_I = Function(FORM_names['cal_I'], *elementary_monomials)

            # initialize the Product to be passed to the subtraction
            subtraction_initializer = Product(monomials, pole_part_initializer, symbolic_cal_I, copy=False)

            # call the subtraction routines
            # Integrate by parts until the `ibp_power_goal` is reached,
            # then do the original subtraction (`integrate_pole_part`).
            at_most_log_poles = integrate_by_parts(subtraction_initializer, ibp_power_goal, *integration_variable_indices)
            subtracted = []
            for item in at_most_log_poles:
                subtracted.extend(  integrate_pole_part(item, *integration_variable_indices)  )

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
                expansion_orders = required_orders + highest_poles_current_term
                regular_expanded = expand_Taylor(regular, regulator_indices, expansion_orders)

                if i == 0: # first iteration; ``highest_poles_current_sector`` not yet set
                    highest_poles_current_sector = highest_poles_current_term
                else:
                    highest_poles_current_sector = np.maximum(highest_poles_current_sector, highest_poles_current_term)

                integrand_summands.append( Product(singular_expanded,regular_expanded,copy=False) )

            integrand = Sum(*integrand_summands, copy=False)

            # update the global lowest
            lowest_orders = np.minimum(lowest_orders, -highest_poles_current_sector)

            # initialize the CFunctions for FORM
            cal_I_derivative_functions = []
            contourdef_Jacobian_derivative_functions = []
            deformed_integration_variable_derivative_functions = []
            decomposed_polynomial_derivatives = []
            other_functions = []

            # compute the required derivatives
            cal_I_derivatives = {}
            other_derivatives = {}
            decomposed_derivatives = {}
            contourdef_Jacobian_derivatives = {}
            deformed_integration_variable_derivatives = {}

            # python dictionaries are unordered but some insertions depend on others --> need an ordering
            ordered_cal_I_derivative_names = []
            ordered_other_derivative_names = []
            ordered_decomposed_derivative_names = []
            ordered_contourdef_Jacobian_derivative_names = []
            ordered_deformed_integration_variable_derivative_names = []

            def update_derivatives(basename, derivative_tracker, full_expression, derivatives=other_derivatives, ordered_derivative_names=ordered_other_derivative_names, functions=other_functions):
                derivatives[basename] = full_expression # include undifferentiated expression
                ordered_derivative_names.append(basename)
                functions.append(basename) # define the symbol as CFunction in FORM
                for multiindex,expression in derivative_tracker.compute_derivatives(full_expression).items():
                    name = _derivative_muliindex_to_name(basename, multiindex)
                    derivatives[name] = expression
                    ordered_derivative_names.append(name)
                    functions.append(name) # define the symbol as CFunction in FORM

            #  - for cal_I
            update_derivatives(basename=FORM_names['cal_I'], derivative_tracker=symbolic_cal_I, full_expression=cal_I,
                               derivatives=cal_I_derivatives, ordered_derivative_names=ordered_cal_I_derivative_names,
                               functions=cal_I_derivative_functions)

            #  - for the `remainder_expression`
            update_derivatives(basename=FORM_names['remainder_expression'],
                               derivative_tracker=symbolic_remainder_expression,
                               full_expression=this_primary_sector_remainder_expression)

            #  - for the additional factors
            if contour_deformation_polynomial is not None:
                update_derivatives(
                    FORM_names['additional_deformation_factor'], # basename
                    symbolic_additional_deformation_factor, # derivative tracker
                    additional_deformation_factor # full_expression
                )
                for k,(symbolic_factor,factor) in \
                enumerate(zip(symbolic_deformation_factors,deformation_factors)):
                    update_derivatives(
                        FORM_names['additional_deformation_factor'] + str(integration_variables[k]), # basename
                        symbolic_factor, # derivative tracker
                        factor, # full_expression
                        deformed_integration_variable_derivatives, # derivatives
                        ordered_deformed_integration_variable_derivative_names, # ordered_derivative_names
                        deformed_integration_variable_derivative_functions # functions
                    )

            #  - for the `other_polynomials`
            for prod, exponentiated_function, basename in zip(sector.other, symbolic_other_polynomials, names_other_polynomials):
                _, expression = prod.factors
                expression.exponent = 1 # exponent is already part of the `tracker`
                update_derivatives(
                    basename=basename, # name as defined in `polynomial_names` or dummy name
                    derivative_tracker=exponentiated_function.base,
                    full_expression=expression
                )

            #  - for the `polynomials_to_decompose`
            for prod, exponentiated_function, basename in zip(sector.cast , symbolic_polynomials_to_decompose, names_polynomials_to_decompose):
                _, expression = prod.factors
                expression.exponent = 1 # exponent is already part of the `tracker`
                update_derivatives(
                    basename=basename, # name as defined in `polynomial_names` or dummy name
                    derivative_tracker=exponentiated_function.base,
                    full_expression=expression,
                    derivatives=decomposed_derivatives,
                    ordered_derivative_names=ordered_decomposed_derivative_names,
                    functions=decomposed_polynomial_derivatives
                )

            if contour_deformation_polynomial is not None:
            #  - for the contour deformation Jacobian
                update_derivatives(
                    FORM_names['contourdef_Jacobian'], # basename
                    symbolic_contourdef_Jacobian, # derivative tracker
                    contourdef_Jacobian_determinant, # full_expression
                    contourdef_Jacobian_derivatives, # derivatives
                    ordered_contourdef_Jacobian_derivative_names, # ordered_derivative_names
                    contourdef_Jacobian_derivative_functions # functions
                )

            #  - for the deformed integration variables
                for undeformed_name,deformed_variable,derivative_tracker in zip(integration_variables,deformed_integration_parameters,symbolic_deformed_variables):
                    update_derivatives(
                        FORM_names['deformed_variable'] + str(undeformed_name), # basename
                        derivative_tracker,
                        deformed_variable, # full_expression
                        deformed_integration_variable_derivatives, # derivatives
                        ordered_deformed_integration_variable_derivative_names, # ordered_derivative_names
                        deformed_integration_variable_derivative_functions # functions
                    )

            #  - for the contour deformation polynomial
                full_expression = sector.cast[contour_deformation_polynomial_index].factors[1]
                full_expression.exponent = 1 # exponent is already part of the `tracker`
                update_derivatives(
                    str(contour_deformation_polynomial), # basename
                    symbolic_contour_deformation_polynomial, # derivative tracker
                    full_expression, # full_expression
                    decomposed_derivatives, # derivatives
                    ordered_decomposed_derivative_names, # ordered_derivative_names
                    decomposed_polynomial_derivatives # functions
                )

            # determine which derivatives of the user input ``functions`` are needed and
            # generate the corresponding c++ "function_declarations"
            for call in function_calls:
                number_of_arguments = call.number_of_arguments
                derivative_symbols = call.derivative_symbols
                functions.update(derivative_symbols)
                for derivative_symbol in derivative_symbols:
                    function_declarations.add( _make_CXX_function_declaration(derivative_symbol, number_of_arguments) )
            other_functions.extend(functions)

            # remove repetitions in `decomposed_polynomial_derivatives`
            decomposed_polynomial_derivatives = set(decomposed_polynomial_derivatives)

            # generate the function definitions for the insertion in FORM
            if contour_deformation_polynomial is not None:
                FORM_deformed_integration_variable_definitions = ''.join(
                    _make_FORM_function_definition(
                        name, deformed_integration_variable_derivatives[name], integration_variables, limit=10**6
                    )
                    for name in ordered_deformed_integration_variable_derivative_names
                )
                FORM_contourdef_Jacobian_derivative_definitions = ''.join(
                    _make_FORM_function_definition(
                        name, contourdef_Jacobian_derivatives[name], integration_variables, limit=10**6
                    )
                    for name in ordered_contourdef_Jacobian_derivative_names
                )
            FORM_cal_I_definitions = ''.join(
                _make_FORM_function_definition(name, cal_I_derivatives[name], symbols_other_polynomials, limit=10**6)
                for name in ordered_cal_I_derivative_names
            )
            FORM_other_definitions = ''.join(
                _make_FORM_function_definition(name, other_derivatives[name], symbols_remainder_expression, limit=10**6)
                for name in ordered_other_derivative_names
            )
            FORM_decomposed_definitions = ''.join(
                _make_FORM_function_definition(name, decomposed_derivatives[name], symbols_remainder_expression, limit=10**6)
                for name in ordered_decomposed_derivative_names
            )

            # generate list over all occuring orders in the regulators
            regulator_powers = list( rangecomb(np.zeros_like(required_orders), required_orders + highest_poles_current_sector) )
            number_of_orders = len(regulator_powers)

            # generate the definitions of the FORM preprocessor variables "shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex'"
            regulator_powers = _make_FORM_shifted_orders(regulator_powers)

            # parse template file "sector.h"
            template_replacements['functions'] = _make_FORM_list(other_functions)
            template_replacements['cal_I_derivatives'] = _make_FORM_list(cal_I_derivative_functions)
            template_replacements['decomposed_polynomial_derivatives'] = _make_FORM_list(decomposed_polynomial_derivatives)
            template_replacements['insert_cal_I_procedure'] = FORM_cal_I_definitions
            template_replacements['insert_other_procedure'] = FORM_other_definitions
            template_replacements['insert_decomposed_procedure'] = FORM_decomposed_definitions
            template_replacements['integrand_definition_procedure'] = _make_FORM_function_definition(internal_prefix+'sDUMMYIntegrand', integrand, args=None, limit=10**6)
            template_replacements['sector_container_initializer'] = _make_CXX_Series_initialization(regulators, -highest_poles_current_sector,
                                                                                                    required_orders, sector_index,
                                                                                                    contour_deformation_polynomial is not None)
            template_replacements['highest_regulator_poles'] = _make_FORM_list(highest_poles_current_sector)
            template_replacements['regulator_powers'] = regulator_powers
            template_replacements['number_of_orders'] = number_of_orders
            parse_template_file(os.path.join(template_sources, 'codegen', 'sector.h'), # source
                                os.path.join(name,             'codegen', 'sector%i.h' % sector_index), # dest
                                template_replacements)

            if contour_deformation_polynomial is not None:
                # parse template file "contour_deformation.h"
                template_replacements['contourdef_Jacobian_derivative_functions'] = _make_FORM_list(contourdef_Jacobian_derivative_functions)
                template_replacements['deformed_integration_variable_derivative_functions'] = _make_FORM_list(deformed_integration_variable_derivative_functions)
                template_replacements['contour_deformation_polynomial'] = contour_deformation_polynomial
                template_replacements['positive_polynomials'] = _make_FORM_list(positive_polynomials)
                template_replacements['insert_deformed_integration_variables_procedure'] = FORM_deformed_integration_variable_definitions
                template_replacements['insert_contourdef_Jacobian_derivatives_procedure'] = FORM_contourdef_Jacobian_derivative_definitions
                template_replacements['deformation_parameters'] = _make_FORM_list(deformation_parameters)
                parse_template_file(os.path.join(template_sources, 'codegen', 'contour_deformation.h'), # source
                                    os.path.join(name,             'codegen', 'contour_deformation_sector%i.h' % sector_index), # dest
                                    template_replacements)

    # expand the `prefactor` to the required orders
    print('expanding the prefactor')
    required_prefactor_orders = requested_orders - lowest_orders
    expanded_prefactor = expand_sympy(prefactor, regulators, required_prefactor_orders)

    # pack the `prefactor` into c++ function that returns a nested `Series`
    # and takes the `real_parameters` and the `complex_parameters`
    prefactor_type = 'secdecutil::Series<' * len(regulators) + 'integrand_return_t' + '>' * len(regulators)
    prefactor_function_body = _make_prefactor_function(expanded_prefactor, real_parameters, complex_parameters)

    # define the return type of "make_integrands"
    make_integrands_return_t = 'std::vector<' + 'secdecutil::Series<' * len(regulators) + \
                               'secdecutil::IntegrandContainer<integrand_return_t, real_t const * const' + \
                               '>' * (len(regulators) + 2)

    # parse the template files "integrands.cpp", "name.hpp", "pole_structures.cpp", "prefactor.cpp", and "functions.hpp"
    template_replacements['function_declarations'] = '\n'.join(function_declarations)
    template_replacements['make_integrands_return_t'] = make_integrands_return_t
    template_replacements['prefactor_type'] = prefactor_type
    template_replacements['prefactor_function_body'] = prefactor_function_body
    template_replacements['number_of_sectors'] = sector_index
    template_replacements['lowest_orders'] = _make_FORM_list(lowest_orders)
    template_replacements['highest_orders'] = _make_FORM_list(required_orders)
    template_replacements['lowest_prefactor_orders'] = _make_FORM_list(-highest_prefactor_pole_orders)
    template_replacements['highest_prefactor_orders'] = _make_FORM_list(required_prefactor_orders)
    template_replacements['sector_includes'] = ''.join( '#include "sector_%i.hpp"\n' % i for i in range(1,sector_index+1) )
    template_replacements['sectors_initializer'] = ','.join( 'integrand_of_sector_%i' % i for i in range(1,sector_index+1) )
    template_replacements['pole_structures_initializer'] = str(pole_structures).replace(' ','').replace('[','{').replace(']','}')
    parse_template_file(os.path.join(template_sources, 'name.hpp'), # source
                        os.path.join(name,            name + '.hpp'), # dest
                        template_replacements)
    for filename in ['integrands.cpp', 'prefactor.cpp', 'pole_structures.cpp', 'functions.hpp']:
        parse_template_file(os.path.join(template_sources, 'src', filename),
                            os.path.join(name,             'src', filename),
                            template_replacements)

    print('"' + name + '" done')

    # print message how to implement the dummy functions if applicable
    have_dummy_functions = True if functions else False
    if have_dummy_functions:
        print(
                 "Declarations of the `functions` and their required derivatives are provided\n" + \
                 "in the file 'src/functions.hpp'. Please refer to that file for further\n" + \
                 "instructions."
             )

    # return the replacements in the template files
    return template_replacements
