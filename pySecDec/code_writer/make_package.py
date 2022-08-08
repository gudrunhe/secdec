"""
This module implements the main program - the function
:func:`.make_package`.

"""

from ..metadata import version, git_id
from ..misc import sympify_symbols, rangecomb, make_cpp_list, chunks
from ..algebra import _Expression, Expression, Polynomial, \
                      ExponentiatedPolynomial, Pow, Product, \
                      ProductRule, Function, Sum, sympify_expression
from .. import decomposition
from ..matrix_sort import iterative_sort, Pak_sort, light_Pak_sort
from ..subtraction import integrate_pole_part, integrate_by_parts, pole_structure as compute_pole_structure
from ..expansion import expand_singular, expand_Taylor, expand_sympy, OrderError
from ..misc import lowest_order, parallel_det, det
from .template_parser import validate_pylink_qmc_transforms, generate_pylink_qmc_macro_dict, parse_template_file, parse_template_tree
from itertools import chain, repeat
from multiprocessing import Pool
from time import strftime
from re import match
from .. import formset
from collections import namedtuple
import inspect
import numpy as np
import sympy as sp
import sys, os
import pySecDecContrib

# The only public object this module provides is the function `make_package`.
# The module is organized in multiple sections dedicated to specific tasks
# to be addressed while writing the c++ package.
# To find the main function `make_package`, it is easiest to full-text search
# for "def make_package(".

_sympy_zero = sympify_expression(0)
_sympy_one = sympify_expression(1)

# sympy symbols are no longer callable starting from version 1.3
_to_function = lambda x: sp.Function(str(x))

# define the internal names to be used in FORM
internal_prefix = 'SecDecInternal'
FORM_names = dict(
    cal_I=internal_prefix+'CalI',
    real_part=internal_prefix+'RealPart',
    cast_polynomial=internal_prefix+'DecoPoly',
    other_polynomial=internal_prefix+'OtherPoly',
    remainder_expression=internal_prefix+'Remainder',
    error_token=internal_prefix+'ErrorToken',
    contourdef_transform=internal_prefix+'Deformation',
    contourdef_Jacobian=internal_prefix+'CondefJac',
    deformed_variable=internal_prefix+'Deformed',
    deformation_parameter_i=internal_prefix+'Lambda%i',
    additional_deformation_factor=internal_prefix+'CondefFac'
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
            if type(expression) is ExponentiatedPolynomial:
                expression.exponent = sympify_expression(str(expression.exponent))
                expression.coeffs = np.array([ sympify_expression(str(coeff)) for coeff in expression.coeffs ])
            elif type(expression) is Polynomial:
                expression = ExponentiatedPolynomial(expression.expolist,
                                                     np.array([ sympify_expression(str(coeff)) for coeff in expression.coeffs ]),
                                                     _sympy_one, # exponent
                                                     polysymbols, copy=False)
            else:
                expression = sympify_expression(str(expression))
                if expression.is_Pow:
                    assert len(expression.args) == 2
                    expression_base = Polynomial.from_expression(expression.args[0], polysymbols)
                    expression_exponent = sympify_expression(str(expression.args[1]))
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
            expression, functions = Expression(expression, polysymbols, follow_functions=True)
        else:
            raise RuntimeError('`_parse_expressions` is only implemented for `target_type`s in %s, not for %s. (raised while parsing `%s`)' \
                                % (set([_Expression, ExponentiatedPolynomial]), target_type, name_of_make_argument_being_parsed))
        parsed.append(expression if target_type is ExponentiatedPolynomial else (expression,functions))
    return parsed

def _validate(name, allow_underscore=False):
    '''
    Check validity of `name` for usage in FORM,
    c++, and the file system. If `allow_underscore`
    is True the underscores are allowed.
    Restrictions are as follows:
     o FORM: no underscores
     o FORM, c++: the first character must be
       in [A-Z, a-z]; i.e. no leading numbers
     o all: no special characters; i.e. only
       characters from the set [A-Z, a-z, 0-9]
     o `name` must not begin with the internal
       prefix

    '''
    if match(r'^[A-Z,a-z,_]+[A-Z,a-z,0-9,_]*$', name) is None:
        raise NameError('"%s" cannot be used as symbol' % name)

    if not allow_underscore and '_' in name:
        raise NameError('"%s" cannot contain an underscore character "_"' % name)

    # disallowed keywords
    banwords = set([
                   # c++ keywords
                   'std','alignas','alignof','and','and_eq','asm','atomic_cancel','atomic_commit',
                   'atomic_noexcept','auto','bitand','bitor','bool','break','case','catch','char',
                   'char16_t','char32_t','class','compl','concept','const','constexpr','const_cast',
                   'continue','decltype','default','delete','do','double','dynamic_cast','else',
                   'enum','explicit','export','extern','false','float','for','friend','goto','if',
                   'import','inline','int','long','module','mutable','namespace','new','noexcept',
                   'not','not_eq','nullptr','operator','or','or_eq','private','protected','public',
                   'register','reinterpret_cast','requires','return','short','signed',
                   'sizeof','static','static_assert','static_cast','struct','switch','synchronized',
                   'template','this','thread_local','throw','true','try','typedef','typeid','typename',
                   'union','unsigned','using','virtual','void','volatile','wchar_t','while','xor','xor_eq',
                   '_Pragma',
                   # c keywords
                   'auto','break','case','char','const','continue','default','do','double','else','enum',
                   'extern','float','for','goto','if','inline','int','long','register','restrict','return',
                   'short','signed','sizeof','static','struct','switch','typedef','union','unsigned','void',
                   'volatile','while','_Alignas','_Alignof','_Atomic','_Bool','_Complex','_Generic',
                   '_Imaginary','_Noreturn','_Static_assert','_Thread_local',
                   'alignas','alignof','bool','complex','imaginary','noreturn','static_assert',
                   'thread_local',
                   # math functions and symbols
                   'log','exp','real','imag','sqrt','I','i_',
                   # Cuba integrator library
                   'cubareal','Vegas','llVegas','Suave','llSuave','Cuhre','llCuhre','Divonne','llDivonne',
                   'integrand_t','peakfinder_t',
               ])
    for banword in banwords:
        if name == banword:
            raise NameError('"%s" cannot be used as symbol' % name)

    # disallowed prefixes
    banstarts = set([internal_prefix,'cuba','atomic','_'])
    for banstart in banstarts:
        if name.lower().startswith(banstart.lower()):
            raise NameError( '"%s" cannot be used as symbol (must not begin with "%s")' % (name,banstart) )

def _convert_input(name, integration_variables, ibp_power_goal, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   real_parameters, complex_parameters, form_optimization_level,
                   form_work_space, form_memory_use, form_threads, form_insertion_depth,
                   contour_deformation_polynomial, positive_polynomials, decomposition_method, pylink_qmc_transforms):
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
        _validate( str(symbol) , allow_underscore=False )
    _validate( str(name) , allow_underscore=True )

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

    # the exponents must be polynomials in the regulators
    for poly in polynomials_to_decompose + other_polynomials:
        try:
            Polynomial.from_expression(poly.exponent, regulators)
        except Exception as error:
            error.args = tuple(['The exponents of the `polynomials_to_decompose` and the `other_polynomials` must be polynomials in the regulators. Error while checking: "%s"' % poly] + [arg for arg in error.args])
            raise
        exponent = Polynomial.from_expression(poly.exponent, integration_variables)
        assert (exponent.expolist == 0).all(), 'The exponents of the `polynomials_to_decompose` and the `other_polynomials` must not depend on the `integration_variables`. Error while checking: "%s"' % poly

    # the `polynomials_to_decompose` must not depend on the `regulators`
    regulator_indices = [i + len(integration_variables) for i in range(len(regulators))]
    for poly in polynomials_to_decompose:
        assert (poly.expolist[:,regulator_indices] == 0).all(), 'The `polynomials_to_decompose` must not depend on the `regulators`. Error while checking: "%s"' % poly

    # convert ``prefactor`` to sympy expression
    prefactor = sympify_expression(prefactor)

    # convert ``requested_orders`` to numpy array
    requested_orders = np.array(requested_orders)

    # convert ``form_threads`` to integer and check positivity
    assert form_threads == int(form_threads), '`form_threads` must be an integer.'
    assert form_threads >= 1, '`form_threads` must not be positive.'

    # compute form.set parameters
    form_setup = formset.Setup((4, 2, 0))
    form_setup.workspace = formset.parse_number(form_work_space)
    form_setup.threads = form_threads
    if form_memory_use is not None:
        requested_memory_use = formset.parse_number(form_memory_use)
        form_setup = form_setup.scale(requested_memory_use, lowest_scale=1.0)
        obtained_memory_use = form_setup.calc()
        if obtained_memory_use > requested_memory_use:
            print( 'warning: FORM memory usage will be limited to ~%s (not %s)' % (
                formset.round_human_readable(obtained_memory_use, True, True),
                form_memory_use))

    # convert ``form_insertion_depth`` to integer and check nonnegativity
    assert form_insertion_depth == int(form_insertion_depth), '`form_insertion_depth` must be an integer.'
    assert form_insertion_depth >= 0, '`form_insertion_depth` must not be negative.'

    # convert ``ibp_power_goal`` to a dictionary between integration variables and the power goal to keep track after primary decomposition
    # layout as "d = {x:-1, y:0}"
    if np.iterable(ibp_power_goal):
        if not isinstance(ibp_power_goal,list):
            ibp_power_goal = list(ibp_power_goal)
        assert len(ibp_power_goal) == len(integration_variables), 'The number of `ibp_power_goal` (%i) must equal the number of integration_variables (%i).' % (len(ibp_power_goal), len(integration_variables))
    else:
        ibp_power_goal = repeat(ibp_power_goal)
    ibp_power_goal = {iv: power_goal for iv, power_goal in zip(integration_variables, ibp_power_goal)}

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

    pylink_qmc_transforms = validate_pylink_qmc_transforms(pylink_qmc_transforms)

    return (name, integration_variables, ibp_power_goal, regulators,
            requested_orders, polynomials_to_decompose, polynomial_names,
            other_polynomials, prefactor, remainder_expression, functions, function_calls,
            real_parameters, complex_parameters, form_optimization_level,
            form_setup, form_insertion_depth, contour_deformation_polynomial, positive_polynomials,
            decomposition_method, pylink_qmc_transforms, symbols_polynomials_to_decompose, symbols_other_polynomials,
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
                                        geometric_ku=        dict(
                                                                     primary=decomposition.iterative.primary_decomposition,
                                                                     secondary=lambda sector, indices: decomposition.geometric.geometric_decomposition_ku(sector, indices, normaliz, workdir)
                                                                 ),
                                        geometric_no_primary=dict(
                                                                     primary=lambda sector, indices: [sector], # no primary decomposition
                                                                     secondary=lambda sector, indices: decomposition.geometric.geometric_decomposition_ku(sector, indices, normaliz, workdir)
                                                                 ),
                                        geometric_infinity_no_primary=dict(
                                                                     primary=lambda sector, indices: [sector], # no primary decomposition
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
                            form_setup, form_insertion_depth, requested_orders,
                            contour_deformation_polynomial, nested_series_type,
                            enforce_complex):
    '''
    Create the `target_directory` (given by `name`) and return the two
    optional arguments passed to :func:`parse_template_tree`.

    '''
    # initialize template replacements
    template_replacements = dict(
                                     name = name,
                                     name_as_char_pack = ','.join([("'" + char + "'") for char in name]),
                                     number_of_real_parameters = len(real_parameters),
                                     real_parameters = _make_FORM_list(real_parameters),
                                     names_of_real_parameters = make_cpp_list(real_parameters),
                                     number_of_complex_parameters = len(complex_parameters),
                                     complex_parameters = _make_FORM_list(complex_parameters),
                                     names_of_complex_parameters = make_cpp_list(complex_parameters),
                                     have_complex_parameters = len(complex_parameters) > 0,
                                     number_of_regulators = len(regulators),
                                     regulators = _make_FORM_list(regulators),
                                     names_of_regulators = make_cpp_list(regulators),
                                     polynomial_names = _make_FORM_list(regulators),
                                     form_optimization_level = form_optimization_level,
                                     form_work_space = form_setup.workspace,
                                     form_small_size = form_setup.smallsize,
                                     form_large_size = form_setup.largesize,
                                     form_terms_in_small = form_setup.termsinsmall,
                                     form_scratch_size = form_setup.scratchsize,
                                     form_sort_io_size = form_setup.sortiosize,
                                     form_threads = form_setup.threads,
                                     form_insertion_depth = form_insertion_depth,
                                     contour_deformation = int(contour_deformation_polynomial is not None),
                                     requested_orders = _make_FORM_list(requested_orders),
                                     nested_series_type = nested_series_type,
                                     pySecDec_version = version,
                                     python_version = sys.version,
                                     pySecDec_git_id = git_id,
                                     contrib_dirname = pySecDecContrib.dirname,
                                     date_time = strftime("%a %d %b %Y %H:%M"),
                                     enforce_complex_return_type=int(bool(enforce_complex)) # make sure that this is either ``0`` or ``1``
                                )

    # configure template parser
    file_renamings = {
                          # replace "name" by the name of the integral
                          'name' : name,
                          'integrate_name.cpp' : 'integrate_' + name + '.cpp',
                          'cuda_integrate_name.cpp' : 'cuda_integrate_' + name + '.cpp',

                          # the files below are specific for each sector --> do not parse globally
                          'contour_deformation.h' : None,
                          'sector.h' : None,
                          'sector.d' : None,

                          # the files below can only be written after the decomposition is completed
                          'integrands.cpp' : None,
                          'name.hpp' : None,
                          'prefactor.cpp' : None,
                          'pole_structures.cpp' : None,
                          'functions.hpp' : None,
                          'pylink.cpp' : None
                     }
    # needs `number_of_sectors` --> can only be written after the decomposition is completed
    # required qmc template instantiations still also need to be determined
    file_renamings['pylink.cpp'] = None
    file_renamings['qmc_template_instantiations.cpp'] = None


    # the files below are only relevant for contour deformation --> do not parse if deactivated
    if contour_deformation_polynomial is None:
        for filename in ['contour_deformation.h','write_contour_deformation.frm']:
            file_renamings[filename] = None

    # get path to the directory with the template files (path relative to directory with this file: "./templates/")
    template_sources = os.path.join(os.path.split(os.path.abspath(__file__))[0],'templates','make_package')

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
        expression.clear_cache()

        if len(str_expression) <= limit: # no need to split
            codelines.append( ("  Id %s = %s;\n" % (name+FORM_args_left_hand_side,str_expression)).replace("**","^") )
            return

        if type(expression) is ProductRule:
            expression = expression.to_sum()

        if type(expression) is Sum:
            del str_expression
            name_next_level = internal_prefix+'fDUMMY'+name+'Part'
            codelines.append(
                "  Id %s = %s;\n" % (
                    name+FORM_args_left_hand_side,
                    '+'.join( name_next_level+str(i)+FORM_args_right_hand_side for i in range(len(expression.summands)) )
                )
            )
            for i,summand in enumerate(expression.summands):
                recursion(name_next_level+str(i), summand)
            return

        if type(expression) is Product:
            del str_expression
            name_next_level = internal_prefix+'fDUMMY'+name+'Part'
            codelines.append(
                "  Id %s = %s;\n" % (
                    name+FORM_args_left_hand_side,
                    '*'.join( name_next_level+str(i)+FORM_args_right_hand_side for i in range(len(expression.factors)) )
                )
            )
            for i,factor in enumerate(expression.factors):
                recursion(name_next_level+str(i), factor)
            return

        if type(expression) is Polynomial:
            del str_expression
            old_n_coeffs = len(expression.coeffs)
            if old_n_coeffs == 1:
                # we need to break down the coefficient, instead of breaking the polynomial into sum terms
                name_next_level = internal_prefix+'fDUMMY'+name+'Part'
                codelines.append(
                    "  Id %s = %s * (%s);\n" % (
                        name+FORM_args_left_hand_side,
                        name_next_level+str(0)+FORM_args_right_hand_side, Polynomial(expression.expolist, [1], expression.polysymbols, copy=False) )
                    )
                recursion(name_next_level+str(0), expression.coeffs[0])
                return
            else:
                # break the polynomial into roughly 10 equal sized terms
                name_next_level = internal_prefix+'fDUMMY'+name+'Part'
                new_n_coeffs = old_n_coeffs//10
                if new_n_coeffs == 0: new_n_coeffs = 1
                new_n_polys = int(np.ceil(old_n_coeffs/new_n_coeffs))
                codelines.append(
                    "  Id %s = %s;\n" % (
                        name+FORM_args_left_hand_side,
                        '+'.join( name_next_level+str(i)+FORM_args_right_hand_side for i in range(new_n_polys) )
                    )
                )
                for i in range(new_n_polys):
                    recursion(name_next_level+str(i), Polynomial(expression.expolist[i*new_n_coeffs:(i+1)*new_n_coeffs,:], expression.coeffs[i*new_n_coeffs:(i+1)*new_n_coeffs], expression.polysymbols, False) )
                return

        # rescue: print warning and write unsplit expression
        print( 'WARNING: Could not split "%s" (not implemented for %s)' % (name,type(expression)) )
        codelines.append( ("  Id %s = %s;\n" % (name+FORM_args_left_hand_side,str_expression)).replace("**","^") )
        return

    recursion(name, expression)
    return codelines

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

def _make_sector_cpp_files(sector_index, regulator_powers, highest_poles, contour_deformation):
    """
    Produce a Makefile-formatted list of .cpp files that
    write_integrand.frm will produce for a given sector.

    Please keep this synchronized with write_integrand.frm,
    because the logic is duplicated.
    """
    files = sorted(
        ("sector_%d_" % sector_index) + "_".join(str(pol-hi) for pol, hi in zip(p, highest_poles)).replace("-", "n") + ".cpp"
        for p in regulator_powers
    )
    if contour_deformation:
        files += ["contour_deformation_" + f for f in files] + ["optimize_deformation_parameters_" + f for f in files]
    files = ["sector_%d.cpp" % sector_index] + files
    return " \\\n\t".join("src/" + f for f in files)

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
                outstr_body_snippets.append( sp.printing.ccode(coeff.evalf(20)) )
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

class MaxDegreeFunction(Function):
    '''
    Symbolic function where derivatives
    higher than `maxdegrees` are assumed to
    be zero.

    '''
    @staticmethod
    def get_maxdegrees(expression, ignore_subclass, indices=None):
        if indices is None:
            indices = range(expression.number_of_variables)
        if (isinstance(expression, Polynomial) if ignore_subclass else type(expression) is Polynomial):
            maxdegrees = np.zeros_like(expression.expolist[0]) + np.inf
            maxdegrees[indices] = expression.expolist[:,indices].max(axis=0)
            maxdegrees[expression.expolist.min(axis=0) < 0] = np.inf # Taylor series of ``1/x^n`` with ``n > 0`` never truncates
        else:
            maxdegrees = np.array( [np.inf] * expression.number_of_variables )
        return maxdegrees

    def __init__(self, symbol, *arguments, **kwargs):
        maxdegrees = kwargs.pop('maxdegrees', np.array( [np.inf] * arguments[0].number_of_variables ))
        super(MaxDegreeFunction, self).__init__(symbol, *arguments, **kwargs)
        self.maxdegrees = list(maxdegrees)

    def derive(self, index):
        # modified `derive` of the base class
        def super_derive(self, index):
            '''
            Generate the derivative by the parameter indexed `index`.
            The derivative of a function with `symbol` ``f`` by
            some `index` is denoted as ``dfd<index>``.

            :param index:
                integer;
                The index of the paramater to derive by.

            '''
            derivative = self.derivatives[index]
            if derivative is not None:
                return derivative.copy()

            summands = []
            for argindex, arg in enumerate(self.arguments):
                if self.maxdegrees[argindex] <= 0:
                    continue

                if self.differentiated_args[argindex, index] is None:
                    differentiated_arg = arg.derive(index)
                    self.differentiated_args[argindex, index] = differentiated_arg
                else:
                    differentiated_arg = self.differentiated_args[argindex, index]

                # catch ``differentiated_arg == 0``
                if type(differentiated_arg) is Polynomial and (differentiated_arg.coeffs == 0).all():
                    continue

                # generate the multiindex of the requested derivative
                old_multiindex = self.derivative_multiindex
                new_multiindex = list(old_multiindex)
                new_multiindex[argindex] += 1
                new_multiindex = tuple(new_multiindex)
                self.derivative_tracks[new_multiindex] = [argindex, old_multiindex]

                derivative_symbol = 'd' * np.sum(new_multiindex) + self.basename + \
                    ''.join( ('d%i' %argindex) * howmany for argindex,howmany in enumerate(new_multiindex) )

                self.derivative_symbols.add(derivative_symbol)
                new_maxdegrees = self.maxdegrees[:]
                new_maxdegrees[argindex] -= 1
                summands.append(
                                ProductRule(    # chain rule
                                                differentiated_arg,
                                                type(self)(derivative_symbol, *(arg.copy() for arg in self.arguments),
                                                           differentiated_args=self.differentiated_args, copy=False,
                                                           derivative_symbols=self.derivative_symbols, basename=self.basename,
                                                           derivative_multiindex=new_multiindex, derivative_tracks=self.derivative_tracks,
                                                           maxdegrees=new_maxdegrees), #modification: additional argument `new_maxdegrees`
                                                copy=False
                                           )
                               )
            if summands:
                derivative = self.derivatives[index] = Sum(*summands, copy=False)
            else: # return zero packed into a `Polynomial`
                derivative = self.derivatives[index] = Polynomial(np.zeros([1,self.number_of_variables], dtype=int), np.array([0]), self.symbols, copy=False)

            return derivative

        derivative = super_derive(self,index)
        return derivative

    def copy(self):
        copy = super(MaxDegreeFunction, self).copy()
        copy.maxdegrees = list(self.maxdegrees)
        return copy

    def replace(self, index, value, remove=False):
        replacement = super(MaxDegreeFunction, self).replace(index, value, remove)
        replacement.maxdegrees = list(self.maxdegrees)
        return replacement

def _make_environment(original_environment):
    'Prepare the environment for :func:`._process_secondary_sector`.'

    original_environment = original_environment.copy()
    sector_index = original_environment.pop('sector_index')
    secondary_sectors = original_environment.pop('secondary_sectors')

    # remove items that are not needed or cannot be pickled
    original_environment.pop('pool', None)
    original_environment.pop('strategy', None)
    original_environment.pop('original_decomposition_strategies', None)
    original_environment.pop('primary_sectors', None)
    original_environment.pop('primary_sectors_to_consider', None)
    original_environment.pop('primary_decomposition_with_splitting', None)
    original_environment.pop('secondary_decomposition_with_splitting', None)

    for sector in secondary_sectors:
        sector_index += 1
        environment = original_environment.copy()
        environment['sector'] = sector
        environment['sector_index'] = sector_index
        yield environment

def _process_secondary_sector(environment):
    'Function to process the `secondary_sectors` in parallel.'

    # read environment
    sector_index = environment['sector_index']
    sector = environment['sector']
    symbols_other_polynomials = environment['symbols_other_polynomials']
    all_symbols = environment['all_symbols']
    contour_deformation_polynomial = environment['contour_deformation_polynomial']
    str_error_token = environment['str_error_token']
    template_sources = environment['template_sources']
    symbols_polynomials_to_decompose = environment['symbols_polynomials_to_decompose']
    regulator_indices = environment['regulator_indices']
    function_calls = environment['function_calls']
    real_parameters = environment['real_parameters']
    pole_part_initializer = environment['pole_part_initializer']
    positive_polynomials = environment['positive_polynomials']
    elementary_monomials_all_symbols = environment['elementary_monomials_all_symbols']
    integration_variable_indices = environment['integration_variable_indices']
    polynomial_names = environment['polynomial_names']
    name = environment['name']
    lowest_orders = environment['lowest_orders']
    imaginary_unit = environment['imaginary_unit']
    have_dummy_functions = environment['have_dummy_functions']
    decomposition_method = environment['decomposition_method']
    normaliz_executable = environment['normaliz_executable']
    use_iterative_sort = environment['use_iterative_sort']
    use_light_Pak = environment['use_light_Pak']
    use_dreadnaut = environment['use_dreadnaut']
    use_Pak = environment['use_Pak']
    complex_parameters = environment['complex_parameters']
    expolist = environment['expolist']
    use_symmetries = environment['use_symmetries']
    sector_index = environment['sector_index']
    integration_variables = environment['integration_variables']
    required_orders = environment['required_orders']
    file_renamings = environment['file_renamings']
    regulators = environment['regulators']
    all_integration_variables = environment['all_integration_variables']
    pole_structures = environment['pole_structures']
    ibp_power_goal_this_primary_sector = environment['ibp_power_goal_this_primary_sector']
    nested_series_type = environment['nested_series_type']
    form_insertion_depth = environment['form_insertion_depth']
    reversed_polynomial_names = environment['reversed_polynomial_names']
    one = environment['one']
    this_primary_sector_remainder_expression = environment['this_primary_sector_remainder_expression']
    transformations = environment['transformations']
    primary_sector = environment['primary_sector']
    remainder_expression_is_trivial = environment['remainder_expression_is_trivial']
    highest_prefactor_pole_orders = environment['highest_prefactor_pole_orders']
    remainder_expression = environment['remainder_expression']
    other_polynomials = environment['other_polynomials']
    polynomial_zero = environment['polynomial_zero']
    function_declarations = environment['function_declarations']
    names_other_polynomials = environment['names_other_polynomials']
    primary_sector_index = environment['primary_sector_index']
    form_optimization_level = environment['form_optimization_level']
    initial_sector = environment['initial_sector']
    polynomial_one = environment['polynomial_one']
    functions = environment['functions']
    requested_orders = environment['requested_orders']
    str_replaced_remainder_expression = environment['str_replaced_remainder_expression']
    split = environment['split']
    elementary_monomials = environment['elementary_monomials']
    symbols_remainder_expression = environment['symbols_remainder_expression']
    symbol = environment['symbol']
    enforce_complex = environment['enforce_complex']
    polynomials_to_decompose = environment['polynomials_to_decompose']
    error_token = environment['error_token']
    template_replacements = environment['template_replacements']
    prefactor = environment['prefactor']
    if contour_deformation_polynomial is not None:
        contourdef_Jacobian_determinant = environment['contourdef_Jacobian_determinant']
        contourdef_Jacobian = environment['contourdef_Jacobian']
        symbolic_deformed_variables = environment['symbolic_deformed_variables']
        symbolic_deformation_factors = environment['symbolic_deformation_factors']
        symbolic_deformation_factor_names = environment['symbolic_deformation_factor_names']
        contour_deformation_polynomial_index = environment['contour_deformation_polynomial_index']
        maxdegrees_contour_deformation_polynomial = environment['maxdegrees_contour_deformation_polynomial']
        symbolic_deformed_variable_names = environment['symbolic_deformed_variable_names']
        symbolic_contour_deformation_polynomial = environment['symbolic_contour_deformation_polynomial']
        str_contour_deformation_polynomial = environment['str_contour_deformation_polynomial']
        deformation_parameters = environment['deformation_parameters']
        deformed_integration_parameters = environment['deformed_integration_parameters']
        deformation_factors = environment['deformation_factors']

    print('writing FORM files for sector', sector_index)

    def parse_exponents(sector, symbols_polynomials_to_decompose, symbols_other_polynomials):
        #  - in ``sector.cast``
        for product in sector.cast:
            mono, poly = product.factors
            mono.exponent = Polynomial.from_expression(mono.exponent, symbols_polynomials_to_decompose)
            poly.exponent = Polynomial.from_expression(poly.exponent, symbols_polynomials_to_decompose)

        #  - in ``sector.other``
        for poly in sector.other:
            try:
                poly.exponent = Polynomial.from_expression(poly.exponent, symbols_other_polynomials)
            except AttributeError:
                pass # not an exponentiated polynomial --> nothing to do

    def parse_coeffs(sector, symbols_polynomials_to_decompose, symbols_other_polynomials):
        #  - in ``sector.cast``
        for product in sector.cast:
            mono, poly = product.factors
            mono.coeffs = np.array([Expression(coeff, symbols_polynomials_to_decompose) for coeff in mono.coeffs])
            poly.coeffs = np.array([Expression(coeff, symbols_polynomials_to_decompose) for coeff in poly.coeffs])

        #  - in ``sector.other``
        for poly in sector.other:
            poly.coeffs = np.array([Expression(coeff, symbols_polynomials_to_decompose) for coeff in poly.coeffs])

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

    # define symbols for the `polynomials_to_decompose` --> shorter expressions and faster in python
    symbolic_polynomials_to_decompose = []
    symbolic_polynomials_to_decompose_all_symbols_undeformed = []
    names_polynomials_to_decompose = []
    for i in range(len(polynomials_to_decompose)):
        try:
            poly_name = str(polynomial_names[i])
        except IndexError:
            poly_name = FORM_names['cast_polynomial'] + str(i)
        symbolic_polynomials_to_decompose_all_symbols_undeformed.append(
                MaxDegreeFunction(
                                  poly_name,
                                  *elementary_monomials_all_symbols,
                                  maxdegrees=MaxDegreeFunction.get_maxdegrees(
                                                                                  sector.cast[i].factors[1],
                                                                                  ignore_subclass=True
                                                                             )
                                  )
        )
        symbolic_polynomials_to_decompose.append(
            Pow(
                MaxDegreeFunction(
                                  poly_name,
                                  *(elementary_monomials if contour_deformation_polynomial is None else symbolic_deformed_variables),
                                  maxdegrees=MaxDegreeFunction.get_maxdegrees(
                                                                                  sector.cast[i].factors[1],
                                                                                  ignore_subclass=True
                                                                             )
                                  ),
                Expression(sector.cast[i].factors[1].exponent, symbols_polynomials_to_decompose),
                copy = False
            )
        )
        names_polynomials_to_decompose.append(poly_name)

    # convert all coefficients to pySecDec expressions
    parse_coeffs(sector, symbols_polynomials_to_decompose+polynomial_names, symbols_other_polynomials+polynomial_names)

    # remove `polynomial_names` - keep polynomial part symbolic as dummy function:
    #  - from `other_polynomials`
    replacements = []
    for k in range(1, len(reversed_polynomial_names) + 1):
        replacement = symbolic_polynomials_to_decompose_all_symbols_undeformed[len(polynomial_names)-k]
        for i in range(k):
            replacement = replacement.replace(-1, error_token, remove=True)
        replacements.append(replacement)
    for replacement in replacements:
        for i, poly in enumerate(sector.other):
            poly = poly.replace(-1, replacement, remove=True)
            for j, coeff in enumerate(poly.coeffs):
                poly.coeffs[j] = coeff.simplify()
            sector.other[i] = poly

    #  - from `polynomials_to_decompose`
    for i in range(len(polynomials_to_decompose)):
        # no dependence here, just remove the symbols
        prod = sector.cast[i]
        for poly_name in reversed_polynomial_names:
            prod = prod.replace(-1, error_token, remove=True)
        sector.cast[i] = prod

    #  - from `Jacobian`
    Jacobian = sector.Jacobian
    for poly_name in reversed_polynomial_names:
        Jacobian = Jacobian.replace(-1, error_token, remove=True)

    #  - from `this_transformations`
    if not use_symmetries:
        for i in range(len(all_integration_variables)):
            for poly_name in reversed_polynomial_names:
                this_transformations[i] = this_transformations[i].replace(-1, error_token, remove=True)

    # convert all exponents to pySecDec expressions
    parse_exponents(sector, symbols_polynomials_to_decompose, symbols_other_polynomials)

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

    # define symbols for the `other_polynomials --> shorter expressions and faster in python
    symbolic_other_polynomials = [
                                     Pow(
                                         MaxDegreeFunction(
                                                     FORM_names['other_polynomial'] + str(i),
                                                     *(elementary_monomials if contour_deformation_polynomial is None else symbolic_deformed_variables),
                                                     maxdegrees=MaxDegreeFunction.get_maxdegrees(sector.other[i].factors[1], indices=regulator_indices, ignore_subclass=True)
                                         ),
                                         Expression(sector.other[i].factors[1].exponent, symbols_polynomials_to_decompose),
                                         copy = False
                                     )
                                     for i in range(len(other_polynomials))
                                 ]

    # Apply ``this_transformation`` and the contour deformaion (if applicable) to
    # the `remainder_expression` BEFORE taking derivatives.
    # Introduce a symbol for the `remainder_expression` and insert in FORM.
    # Note: ``elementary_monomials[len(integration_variables):]`` are the regulators
    maxdegrees_remainder_expression = MaxDegreeFunction.get_maxdegrees(this_primary_sector_remainder_expression, indices=regulator_indices, ignore_subclass=False)
    if remainder_expression_is_trivial:
        # `remainder_expression` does not depend on the integration variables in that case
        maxdegrees_remainder_expression[:len(integration_variables)] = 0
        symbolic_remainder_expression_arguments = [polynomial_zero] * len(integration_variables) + elementary_monomials[len(integration_variables):]
    else:
        symbolic_remainder_expression_arguments = this_transformations + elementary_monomials[len(integration_variables):]
    symbolic_remainder_expression = MaxDegreeFunction(FORM_names['remainder_expression'], *symbolic_remainder_expression_arguments, maxdegrees=maxdegrees_remainder_expression)

    # initialize the product of monomials for the subtraction
    monomial_factors = list(chain([Jacobian.factors[0]], (prod.factors[0] for prod in sector.cast), (prod.factors[0] for prod in sector.other)))
    monomials = Product(*monomial_factors, copy=False)

    # compute the pole structure
    this_pole_structures = [sp.printing.ccode(ps) for ps in compute_pole_structure(monomials, *integration_variable_indices)]

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
    # don't use symbolic names for constant factors in the cal_I --> derivatives of cal_I much shorter
    def is_constant(poly):
        if isinstance(poly,Polynomial) and not np.any(poly.expolist):
            return True
        return False
    nontrivial_symbolic_polynomials_to_decompose = [poly.factors[1] if is_constant(poly.factors[1]) else symbolic_poly for symbolic_poly,poly in zip(symbolic_polynomials_to_decompose,sector.cast)]
    nontrivial_symbolic_other_polynomials = [poly.factors[1] if is_constant(poly.factors[1]) else symbolic_poly for symbolic_poly,poly in zip(symbolic_other_polynomials,sector.other)]
    nontrivial_remainder_expression = this_primary_sector_remainder_expression if is_constant(this_primary_sector_remainder_expression) else symbolic_remainder_expression
    cal_I = Product(nontrivial_remainder_expression, *chain([Jacobian.factors[1]], nontrivial_symbolic_polynomials_to_decompose, nontrivial_symbolic_other_polynomials), copy=False)

    # multiply Jacobian determinant to `cal_I`
    if contour_deformation_polynomial is not None:
        symbolic_contourdef_Jacobian = Function(FORM_names['contourdef_Jacobian'], *elementary_monomials[:len(integration_variables)])
        symbolic_additional_deformation_factor = Function(FORM_names['additional_deformation_factor'], *elementary_monomials)
        nontrivial_symbolic_contourdef_Jacobian = contourdef_Jacobian_determinant if is_constant(contourdef_Jacobian_determinant) else symbolic_contourdef_Jacobian
        nontrivial_symbolic_additional_deformation_factor = additional_deformation_factor if is_constant(additional_deformation_factor) else symbolic_additional_deformation_factor
        cal_I = Product(nontrivial_symbolic_contourdef_Jacobian, nontrivial_symbolic_additional_deformation_factor, cal_I)

    cal_I.simplify()

    # it is faster to use a dummy function for ``cal_I`` and substitute back in FORM
    symbolic_cal_I = Function(FORM_names['cal_I'], *elementary_monomials)

    # initialize the Product to be passed to the subtraction
    subtraction_initializer = Product(monomials, pole_part_initializer, symbolic_cal_I, copy=False)

    # call the subtraction routines
    # Integrate by parts until the `ibp_power_goal` is reached,
    # then do the original subtraction (`integrate_pole_part`).
    after_ibp = integrate_by_parts(subtraction_initializer, ibp_power_goal_this_primary_sector, integration_variable_indices)
    subtracted = []
    for item in after_ibp:
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

        try:
            singular_expanded = expand_singular(Product(singular, copy=False), regulator_indices, required_orders)
        except OrderError:
            coeffs = np.array([_sympy_zero])
            expolist = np.zeros((1, len(singular.symbols)), dtype=int)
            expolist[:,regulator_indices] = required_orders
            singular_expanded = Polynomial(expolist, coeffs, singular.symbols, copy=False)

        highest_poles_current_term = - singular_expanded.expolist[:,regulator_indices].min(axis=0)
        expansion_orders = required_orders + highest_poles_current_term

        try:
            regular_expanded = expand_Taylor(regular, regulator_indices, expansion_orders)
        except OrderError:
            coeffs = np.array([_sympy_zero])
            expolist = np.zeros((1, len(singular.symbols)), dtype=int)
            expolist[:,regulator_indices] = expansion_orders
            regular_expanded = Polynomial(expolist, coeffs, singular.symbols, copy=False)

        if i == 0: # first iteration; ``highest_poles_current_sector`` not yet set
            highest_poles_current_sector = highest_poles_current_term
        else:
            highest_poles_current_sector = np.maximum(highest_poles_current_sector, highest_poles_current_term)

        integrand_summands.append( Product(singular_expanded,regular_expanded,copy=False) )

    integrand = Sum(*integrand_summands, copy=False)

    # update the `lowest_orders`
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
        _, expression = prod.copy().factors
        expression.exponent = 1 # exponent is already part of the `tracker`
        update_derivatives(
            basename=basename, # name as defined in `polynomial_names` or dummy name
            derivative_tracker=exponentiated_function.base,
            full_expression=expression
        )

    #  - for the `polynomials_to_decompose`
    for i,symbolic_polynomials in enumerate((symbolic_polynomials_to_decompose,symbolic_polynomials_to_decompose_all_symbols_undeformed)):
        for prod, exponentiated_function, basename in zip(sector.cast , symbolic_polynomials, names_polynomials_to_decompose):
            _, expression = prod.copy().factors # take copy to avoid modifying original expression
            expression.exponent = 1 # exponent is already part of the `tracker`
            update_derivatives(
                basename=basename, # name as defined in `polynomial_names` or dummy name
                derivative_tracker=exponentiated_function.base if i == 0 else exponentiated_function,
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
        # deformed variable derivative tracker is copied in ``parallel_det(contourdef_Jacobian, pool)``
        #   --> make sure all required derivatives are generated anyway
        for undeformed_name,deformed_variable in zip(integration_variables,symbolic_deformed_variables):
            for j in range(len(integration_variables)):
                symbolic_contourdef_Jacobian.compute_derivatives(deformed_variable.derive(j))

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
        full_expression = sector.cast[contour_deformation_polynomial_index].copy().factors[1]
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
    ordered_decomposed_derivative_names = set(ordered_decomposed_derivative_names)

    # generate the function definitions for the insertion in FORM
    if contour_deformation_polynomial is not None:
        FORM_vanishing_deformed_integration_variable_calls = ''.join(
                    '  Id %s(' % _derivative_muliindex_to_name(FORM_names['deformed_variable'] + str(outer_var), multiindex) + \
                    ','.join(str(inner_var) + ('?{0,1}' if i == j else '?') for j,inner_var in enumerate(integration_variables)) + \
                    ') = %s;\n' % ('0' if np.any(multiindex) else str(outer_var))
                    if multiindex[i] == 0 else ''
                for i,outer_var in enumerate(integration_variables)
            for multiindex in chain([[0]*len(integration_variables)], symbolic_deformed_variables[i].derivative_tracks.keys())
        )
        FORM_deformed_integration_variable_definitions = list(
            _make_FORM_function_definition(
                name, deformed_integration_variable_derivatives[name], integration_variables, limit=10**6
            )
            for name in ordered_deformed_integration_variable_derivative_names
        )
        FORM_contourdef_Jacobian_derivative_definitions = list(
            _make_FORM_function_definition(
                name, contourdef_Jacobian_derivatives[name], integration_variables, limit=10**6
            )
            for name in ordered_contourdef_Jacobian_derivative_names
        )
    FORM_cal_I_definitions = list(
        _make_FORM_function_definition(name, cal_I_derivatives[name], symbols_other_polynomials, limit=10**6)
        for name in ordered_cal_I_derivative_names
    )
    FORM_other_definitions = list(
        _make_FORM_function_definition(name, other_derivatives[name], symbols_remainder_expression, limit=10**6)
        for name in ordered_other_derivative_names
    )
    FORM_decomposed_definitions = list(
        _make_FORM_function_definition(name, decomposed_derivatives[name], symbols_remainder_expression, limit=10**6)
        for name in ordered_decomposed_derivative_names
    )

    # generate list over all occuring orders in the regulators
    regulator_powers = list( rangecomb(np.zeros_like(required_orders), required_orders + highest_poles_current_sector) )
    number_of_orders = len(regulator_powers)

    # generate the definitions of the FORM preprocessor variables "shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex'"
    sector_cpp_files = _make_sector_cpp_files(sector_index, regulator_powers, highest_poles_current_sector, contour_deformation_polynomial is not None)
    regulator_powers = _make_FORM_shifted_orders(regulator_powers)

    # parse template file "sector.h"
    template_replacements['sector_index'] = sector_index
    template_replacements['functions'] = _make_FORM_list(other_functions)
    template_replacements['cal_I_derivatives'] = _make_FORM_list(cal_I_derivative_functions)
    template_replacements['decomposed_polynomial_derivatives'] = _make_FORM_list(decomposed_polynomial_derivatives)
    template_replacements['insert_cal_I_procedure'] = FORM_cal_I_definitions
    template_replacements['insert_other_procedure'] = FORM_other_definitions
    template_replacements['insert_decomposed_procedure'] = FORM_decomposed_definitions
    template_replacements['integrand_definition_procedure'] = _make_FORM_function_definition(internal_prefix+'sDUMMYIntegrand', integrand, args=None, limit=10**6)
    template_replacements['highest_regulator_poles'] = _make_FORM_list(highest_poles_current_sector)
    template_replacements['required_orders'] = _make_FORM_list(required_orders)
    template_replacements['regulator_powers'] = regulator_powers
    template_replacements['number_of_orders'] = number_of_orders
    template_replacements['sector_cpp_files'] = sector_cpp_files
    template_replacements['sector_hpp_files'] = sector_cpp_files.replace(".cpp", ".hpp")
    template_replacements['sector_codegen_sources'] = \
            "codegen/sector%i.h" % sector_index if contour_deformation_polynomial is None else \
            "codegen/sector%i.h codegen/contour_deformation_sector%i.h" % (sector_index, sector_index)
    parse_template_file(os.path.join(template_sources, 'codegen', 'sector.h'), # source
                        os.path.join(name,             'codegen', 'sector%i.h' % sector_index), # dest
                        template_replacements)
    parse_template_file(os.path.join(template_sources, 'codegen', 'sector.d'), # source
                        os.path.join(name,             'codegen', 'sector%i.d' % sector_index), # dest
                        template_replacements)
    for key in 'functions', 'cal_I_derivatives', 'decomposed_polynomial_derivatives','insert_cal_I_procedure','insert_other_procedure','insert_decomposed_procedure', \
            'integrand_definition_procedure','highest_regulator_poles','required_orders','regulator_powers','number_of_orders', \
            'sector_index', 'sector_cpp_files', 'sector_hpp_files', 'sector_codegen_sources':
        del template_replacements[key]

    if contour_deformation_polynomial is not None:
        # parse template file "contour_deformation.h"
        template_replacements['contourdef_Jacobian_derivative_functions'] = _make_FORM_list(contourdef_Jacobian_derivative_functions)
        template_replacements['deformed_integration_variable_derivative_functions'] = _make_FORM_list(deformed_integration_variable_derivative_functions)
        template_replacements['contour_deformation_polynomial'] = contour_deformation_polynomial
        template_replacements['positive_polynomials'] = _make_FORM_list(positive_polynomials)
        template_replacements['nullify_vanishing_deformed_integration_variable_calls_procedure'] = FORM_vanishing_deformed_integration_variable_calls
        template_replacements['insert_deformed_integration_variables_procedure'] = FORM_deformed_integration_variable_definitions
        template_replacements['insert_contourdef_Jacobian_derivatives_procedure'] = FORM_contourdef_Jacobian_derivative_definitions
        template_replacements['deformation_parameters'] = _make_FORM_list(deformation_parameters)
        parse_template_file(os.path.join(template_sources, 'codegen', 'contour_deformation.h'), # source
                            os.path.join(name,             'codegen', 'contour_deformation_sector%i.h' % sector_index), # dest
                            template_replacements)
        for key in 'contourdef_Jacobian_derivative_functions','deformed_integration_variable_derivative_functions','contour_deformation_polynomial','positive_polynomials', \
                'nullify_vanishing_deformed_integration_variable_calls_procedure','insert_deformed_integration_variables_procedure','insert_contourdef_Jacobian_derivatives_procedure', \
                'deformation_parameters':
            del template_replacements[key]

    return lowest_orders, function_declarations, this_pole_structures

def _reduce_sectors_by_symmetries(sectors, message, indices, use_iterative_sort, use_light_Pak_sort, use_Pak, use_dreadnaut, name):
    '''
    Function that reduces the number of sectors by
    identifying symmetries.

    '''
    print(message + ' before symmetry finding:', len(sectors))
    # find symmetries
    if use_iterative_sort:
        sectors = decomposition.squash_symmetry_redundant_sectors_sort(sectors, iterative_sort, indices)
        print(message + ' after symmetry finding (iterative):', len(sectors))
    if use_light_Pak_sort:
        sectors = decomposition.squash_symmetry_redundant_sectors_sort(sectors, light_Pak_sort, indices)
        print(message + ' after symmetry finding (light Pak):', len(sectors))
    if use_Pak:
        sectors = decomposition.squash_symmetry_redundant_sectors_sort(sectors, Pak_sort, indices)
        print(message + ' after symmetry finding (full Pak):', len(sectors))
    if use_dreadnaut:
        sectors = decomposition.squash_symmetry_redundant_sectors_dreadnaut(sectors, indices, use_dreadnaut, os.path.join(name,'dreadnaut_workdir'))
        print(message + ' after symmetry finding (dreadnaut):', len(sectors))
    return sectors


# ---------------------------------- main function ----------------------------------
def make_package(name, integration_variables, regulators, requested_orders,
                 polynomials_to_decompose, polynomial_names=[], other_polynomials=[],
                 prefactor=1, remainder_expression=1, functions=[], real_parameters=[],
                 complex_parameters=[], form_optimization_level=2, form_work_space='50M',
                 form_memory_use=None, form_threads=2,
                 form_insertion_depth=5, contour_deformation_polynomial=None, positive_polynomials=[],
                 decomposition_method='iterative_no_primary', normaliz_executable='normaliz',
                 enforce_complex=False, split=False, ibp_power_goal=-1, use_iterative_sort=True,
                 use_light_Pak=True, use_dreadnaut=False, use_Pak=True, processes=None, pylink_qmc_transforms=['korobov3x3']):
    r'''
    Decompose, subtract and expand an expression.
    Return it as c++ package.

    .. seealso::
        In order to decompose a loop integral,
        use the function
        :func:`pySecDec.loop_integral.loop_package`.

    .. seealso::
        The generated library is described in
        :ref:`generated_cpp_libs`.

    :param name:
        string;
        The name of the c++ namepace and the output
        directory.

    :param integration_variables:
        iterable of strings or sympy symbols;
        The variables that are to be integrated. The
        intgration region depends on the chosen
        `decomposition_method`.

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
            Only user-defined functions that are provided as
            c++-callable code should be mentioned here.
            Listing basic mathematical functions (e.g. ``log``,
            ``pow``, ``exp``, ``sqrt``, ...) is not required
            and considered an error to avoid name conflicts.

        .. note::
            The power function `pow` and the logarithm `log`
            use the nonstandard continuation with an
            infinitesimal negative imaginary part on the
            negative real axis (e.g. ``log(-1) = -i*pi``).

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as real variables.

    :param complex_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as complex variables.

    :param form_optimization_level:
        integer out of the interval [0,4], optional;
        The optimization level to be used in FORM.
        Default: ``2``.

    :param form_work_space:
        string, optional;
        The FORM WorkSpace. Default: ``'50M'``.

        Setting this to smaller values will reduce FORM memory
        usage (without affecting performance), but each problem
        has some minimum value below which FORM will refuse to
        work: it will fail with error message indicating that
        larger WorkSpace is needed, at which point WorkSpace
        will be adjusted and FORM will be re-run.

    :param form_memory_use:
        string, optional;
        The target FORM memory usage. When specified, `form.set`
        parameters will be adjusted so that FORM uses at most
        approximately this much resident memory.

        The minimum is approximately to 600M + 350M per worker
        thread if ``form_work_space`` is left at ``'50M'``.
        if form_work_space is increased to ``'500M'``, then
        the minimum is 2.5G + 2.5G per worker thread.
        Default: ``None``, meaning use the default FORM values.

    :param form_threads:
        integer, optional;
        Number of threads (T)FORM will use. Default: ``2``.

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

        * 'iterative_no_primary' (default): integration region
          :math:`[0,1]^N`.
        * 'geometric_no_primary': integration region :math:`[0,1]^N`.
        * 'geometric_infinity_no_primary': integration region
          :math:`[0,\infty]^N`.
        * 'iterative': primary decomposition followed by
          integration over :math:`[0,1]^{N-1}`.
        * 'geometric': :math:`x_N` is set to one followed by
          integration over :math:`[0,\infty]^{N-1}`.
        * 'geometric_ku': primary decomposition followed by
          integration over :math:`[0,1]^{N-1}`.

        'iterative', 'geometric', and 'geometric_ku' are only
        valid for loop integrals. An end user should use
        'iterative_no_primary', 'geometric_no_primary', or
        'geometric_infinity_no_primary' here.
        In order to compute loop integrals, please use the
        function :func:`pySecDec.loop_integral.loop_package`.

    :param normaliz_executable:
        string, optional;
        The command to run `normaliz`. `normaliz` is only
        required if `decomposition_method` starts with
        'geometric'.
        Default: 'normaliz'

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
        number or iterable of number, optional;
        The `power_goal` that is forwarded to
        :func:`.integrate_by_parts`.

        This option controls how the subtraction terms are
        generated. Setting it to ``-numpy.inf`` disables
        :func:`.integrate_by_parts`, while ``0`` disables
        :func:`.integrate_pole_part`.

        .. versionadded: 1.4
            A separate power_goal for each of the
            `integration_variables` can be set by passing an
            iterable.

        .. seealso::
            To generate the subtraction terms, this function
            first calls :func:`.integrate_by_parts` for each
            integration variable with the give `ibp_power_goal`.
            Then :func:`.integrate_pole_part` is called.

        Default: ``-1``

    :param use_iterative_sort:
        bool;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_sort`
        with :func:`.iterative_sort` to find sector symmetries.
        Default: ``True``

    :param use_light_Pak:
        bool;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_sort`
        with :func:`.light_Pak_sort` to find sector symmetries.
        Default: ``True``

    :param use_dreadnaut:
        bool or string, optional;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_dreadnaut`
        to find sector symmetries.
        If given a string, interpret that string as the command
        line executable `dreadnaut`. If ``True``, try
        ``$SECDEC_CONTRIB/bin/dreadnaut`` and, if the
        environment variable ``$SECDEC_CONTRIB`` is not set,
        ``dreadnaut``.
        Default: ``False``

    :param use_Pak:
        bool;
        Whether or not to use
        :func:`.squash_symmetry_redundant_sectors_sort`
        with :func:`.Pak_sort` to find sector symmetries.
        Default: ``True``

    :param processes:
        integer or None, optional;
        Parallelize the package generation using at most this many
        processes. If ``None``, use the total number of logical
        CPUs on the system (that is, :func:`os.cpu_count()`),
        or the number of CPUs allocated to the current process
        (``len(os.sched_getaffinity(0))``), on platforms where
        this information is available (i.e. Linux+glibc).
        `New in version 1.3`.
        Default: ``None``

    :param pylink_qmc_transforms:
        list or None, optional;
        Required qmc integral transforms, options are:

        * ``korobov<i>x<j>`` for 1 <= i,j <= 6
        * ``korobov<i>`` for 1 <= i <= 6 (same as ``korobov<i>x<i>``)
        * ``sidi<i>`` for 1 <= i <= 6

        `New in version 1.5`.
        Default: ``['korobov3x3']``
    '''
    print('running "make_package" for "' + name + '"')

    # convert input data types to the data types we need
    name, integration_variables, ibp_power_goal, regulators, \
    requested_orders, polynomials_to_decompose, polynomial_names, \
    other_polynomials, prefactor, remainder_expression, functions, \
    function_calls, real_parameters, complex_parameters, \
    form_optimization_level, form_setup, form_insertion_depth, \
    contour_deformation_polynomial, positive_polynomials, decomposition_method, pylink_qmc_transforms, \
    symbols_polynomials_to_decompose, symbols_other_polynomials, \
    symbols_remainder_expression, all_symbols = \
    _convert_input(name, integration_variables, ibp_power_goal, regulators,
                   requested_orders, polynomials_to_decompose, polynomial_names,
                   other_polynomials, prefactor, remainder_expression, functions,
                   real_parameters, complex_parameters, form_optimization_level,
                   form_work_space, form_memory_use, form_threads, form_insertion_depth,
                   contour_deformation_polynomial, positive_polynomials, decomposition_method, pylink_qmc_transforms)

    # construct the c++ type "nested_series_t"
    # for two regulators, the resulting code should read:
    # "secdecutil::Series<secdecutil::Series<T>>"
    nested_series_type = 'secdecutil::Series<' * len(regulators) + 'T' + '>' * len(regulators)

    # configure the template parser and parse global files
    template_sources, template_replacements, file_renamings = \
        _parse_global_templates(
        name, regulators, polynomial_names,
        real_parameters, complex_parameters, form_optimization_level,
        form_setup, form_insertion_depth, requested_orders,
        contour_deformation_polynomial, nested_series_type,
        enforce_complex
    )

    # get the highest poles from the ``prefactor``
    highest_prefactor_pole_orders = -np.array([lowest_order(prefactor, regulator) for regulator in regulators])

    # compute the required expansion order accounting for the prefactor
    required_orders = requested_orders + highest_prefactor_pole_orders

    # get the decomposition routines
    strategy = get_decomposition_routines(decomposition_method, normaliz_executable, os.path.join(name,'normaliz_workdir'))

    # get dreadnaut command if desired
    if use_dreadnaut:
        if isinstance(use_dreadnaut, str):
            dreadnaut_executable = use_dreadnaut
        else:
            dreadnaut_executable = os.path.join(pySecDecContrib.dirname, 'bin', 'dreadnaut')

    # define the monomials "x0", "x1", "x2", ... to keep track of the transformations
    one = Polynomial([[0]*len(symbols_other_polynomials)], [1], symbols_other_polynomials)
    transformations = []
    for i in range(len(integration_variables)):
        to_append = one.copy()
        to_append.expolist[:,i] = 1
        transformations.append(to_append)

    # define an error token that is multiplied to expressions that should evaluate to zero
    str_error_token = FORM_names['error_token']
    error_token = sp.symbols(str_error_token)

    # make a copy of the `integration_variables` for later reference
    all_integration_variables = list(integration_variables)

    # intialize the c++ declarations of the `functions`
    function_declarations = set()

    have_dummy_functions = True if functions else False

    # initialize the decomposition
    initial_sector = decomposition.Sector(polynomials_to_decompose, other_polynomials + transformations)

    # use iterative method for 1D integrals with Cheng-Wu (avoids problem with 1L tadpole and geometric decomposition)
    if decomposition_method == 'geometric' and len(all_integration_variables) == 1:
        strategy = get_decomposition_routines('iterative', normaliz_executable, os.path.join(name,'normaliz_workdir'))

    # if splitting desired, implement it as additional primary decomposition
    if split:
        # cannot split when using the geometric decomposition method because the integration interval is [0,inf] after the primary decomposition
        if decomposition_method == 'geometric':
            raise ValueError('Cannot have ``split=True`` and ``decomposition_method="geometric"``. You probably want to try ``split=True`` and ``decomposition_method="geometric_ku"``')

        original_decomposition_strategies = strategy

        def primary_decomposition_with_splitting(sector, indices):
            # investigate symmetries before the split
            if use_symmetries:
                primary_sectors = _reduce_sectors_by_symmetries\
                (
                    list(  original_decomposition_strategies['primary'](sector, indices)  ),
                    'number of primary sectors',
                    indices[:-1], # primary decomposition removes one integration variable
                    use_iterative_sort,
                    use_light_Pak,
                    use_Pak,
                    dreadnaut_executable if use_dreadnaut else False,
                    name
                )
            else:
                primary_sectors = original_decomposition_strategies['primary'](sector, indices)
            for output_sector in primary_sectors:
                yield output_sector

        def secondary_decomposition_with_splitting(sector, indices, split_sectors=None):
            if split_sectors is None:
                # split and decompose the `sector`
                split_sectors = decomposition.splitting.split_singular(sector, split, indices)
            for split_sector in split_sectors:
                for decomposed_sector in original_decomposition_strategies['secondary'](split_sector, indices):
                    # check if another split is necessary
                    split_decomposed_sector = list( decomposition.splitting.split_singular(decomposed_sector, split, indices) )
                    if len(split_decomposed_sector) == 1:
                        yield decomposed_sector
                    else:
                        for deeper_sector in secondary_decomposition_with_splitting(decomposed_sector, indices, split_decomposed_sector):
                            yield deeper_sector

        # apply modifications to the decomposition strategy
        strategy = dict(primary=primary_decomposition_with_splitting, secondary=secondary_decomposition_with_splitting)

    # initialize the counter
    sector_index = 0

    # initialize the global lowest orders
    lowest_orders = requested_orders.copy()

    # initialize list of pole structures
    pole_structures = []

    # define the imaginary unit
    imaginary_unit = sympify_expression('I')

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

    # investigate if we can take advantage of sector symmetries
    # We can simplify due to symmetries if the `remainder_expression`
    # does explicitly not depend on any integration variable and if
    # we do not split the integration region.
    remainder_expression_is_trivial = True
    for i in range(len(integration_variables)):
        str_replaced_remainder_expression = str(  remainder_expression.replace(i, error_token)  )
        if str_error_token in str_replaced_remainder_expression:
            remainder_expression_is_trivial = False
            break
    use_symmetries = remainder_expression_is_trivial

    # Check that either the primary decomposition or the `remainder_expression` is trivial.
    # Note that the primary decomposition is specialized for loop integrals.
    if not remainder_expression_is_trivial and decomposition_method not in ('iterative_no_primary', 'geometric_no_primary'):
        raise NotImplementedError('The primary decomposition is only implemented for loop integrals. Please perform the primary decomposition yourself and choose ``decomposition_method="iterative_no_primary"`` or ``decomposition_method="geometric_no_primary"``.')

    if use_symmetries:
        # can investigate sector symmetries
        # we do not need `transformations` in that case --> remove from `initial_sector`
        for i in range(len(integration_variables)):
            initial_sector.other.pop()

    # symmetries are applied elsewhere if we split
    if use_symmetries and not split:
        # run primary decomposition and squash symmetry-equal sectors (using both implemented strategies)
        indices = range(len(integration_variables))
        primary_sectors = list(  strategy['primary'](initial_sector, indices)  )
        if len(primary_sectors) > 1: # no need to look for symmetries if only one sector
            primary_sectors = _reduce_sectors_by_symmetries\
            (
                primary_sectors,
                'number of primary sectors',
                indices[:-1], # primary decomposition removes one integration variable
                use_iterative_sort,
                use_light_Pak,
                use_Pak,
                dreadnaut_executable if use_dreadnaut else False,
                name
            )

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

    # initialize the multiprocessing pool to process the `secondary_sectors` in parallel
    if processes is None:
        try:
            processes = len(os.sched_getaffinity(0))
        except AttributeError:
            processes = os.cpu_count()
    if processes > 1:
        pool = Pool(processes)
    else:
        pool = None

    # try-finally block to make sure that the pool is closed
    try:
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

            # get the ibp power goals for the remaining integration variables
            ibp_power_goal_this_primary_sector = [ibp_power_goal[integration_variables[index]] for index in integration_variable_indices]

            # define "elementary" `_Expression`s such as ``x = Polynomial.from_expression('x', polysymbols)`` for all x
            elementary_monomials = []
            for i, symbol in enumerate(symbols_polynomials_to_decompose):
                expolist = np.zeros([1,len(symbols_polynomials_to_decompose)], dtype=int)
                expolist[:,i] = 1
                elementary_monomials.append( Polynomial(expolist, np.array([1]), symbols_polynomials_to_decompose, copy=False) )
            elementary_monomials_all_symbols = []
            for i, symbol in enumerate(symbols_polynomials_to_decompose):
                expolist = np.zeros([1,len(all_symbols)], dtype=int)
                expolist[:,i] = 1
                elementary_monomials_all_symbols.append( Polynomial(expolist, np.array([1]), all_symbols, copy=False) )

            # we later need ``0`` and ``1`` packed into specific types
            polynomial_zero = Polynomial(np.zeros([1,len(symbols_other_polynomials)], dtype=int), np.array([0]), symbols_other_polynomials, copy=False)
            polynomial_one = Polynomial(np.zeros([1,len(symbols_other_polynomials)], dtype=int), np.array([1]), symbols_other_polynomials, copy=False)
            pole_part_initializer = Pow(polynomial_one, -polynomial_one)

            if contour_deformation_polynomial is not None:
                symbolic_deformed_variable_names = [FORM_names['deformed_variable'] + str(original_varname) for original_varname in integration_variables]
                symbolic_deformation_factor_names = [FORM_names['additional_deformation_factor'] + str(original_varname) for original_varname in integration_variables]

                # Need all first and second derivatives of the `contour_deformation_polynomial`.
                # Since the `contour_deformation_polynomial` is left symbolic they are equal for every subsector after primary decomposition.
                # the maximal degrees in the regulators are unchanged by the decomposition
                maxdegrees_contour_deformation_polynomial = MaxDegreeFunction.get_maxdegrees(primary_sector.cast[contour_deformation_polynomial_index].factors[1], indices=regulator_indices, ignore_subclass=True)
                symbolic_contour_deformation_polynomial = MaxDegreeFunction(str(contour_deformation_polynomial), *elementary_monomials, maxdegrees=maxdegrees_contour_deformation_polynomial)

                # compute the deformation of the integration parameters and its Jacobian matrix (see e.g. section 3.2 in arXiv:1601.03982):
                # ``z_k({x_k}) = x_k * (1 - i * lambda_k * (1-x_k) * Re(dF_dx_k))``, where "dF_dx_k" denotes the derivative of ``F`` by ``x_k``
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
                symbolic_deformed_variables = [
                                                    MaxDegreeFunction(deformed_name,
                                                                      *elementary_monomials[:len(integration_variables)],
                                                                      maxdegrees=maxdegrees_contour_deformation_polynomial)
                                                    for deformed_name in symbolic_deformed_variable_names
                                              ]
                symbolic_deformed_variables.extend( (regulator for regulator in elementary_monomials[len(integration_variables):]) )
                symbolic_deformation_factors = [
                                                    MaxDegreeFunction(deformed_name,
                                                                      *elementary_monomials[:len(integration_variables)],
                                                                      maxdegrees=maxdegrees_contour_deformation_polynomial)
                                                    for deformed_name in symbolic_deformation_factor_names
                                               ]

                # generate the Jacobian determinant
                print('computing Jacobian determinant for primary sector', primary_sector_index)
                contourdef_Jacobian = np.empty((len(integration_variables), len(integration_variables)), dtype=object)
                for i in range(len(integration_variables)):
                    for j in range(len(integration_variables)):
                        contourdef_Jacobian[i,j] = symbolic_deformed_variables[i].derive(j).simplify()
                if pool is not None:
                    contourdef_Jacobian_determinant = parallel_det(contourdef_Jacobian, pool)
                else:
                    contourdef_Jacobian_determinant = det(contourdef_Jacobian)

            # remove `polynomial_names` from the `remainder_expression`
            this_primary_sector_remainder_expression = remainder_expression
            for poly_name in reversed_polynomial_names:
                this_primary_sector_remainder_expression = this_primary_sector_remainder_expression.replace(-1, _to_function(poly_name)(*symbols_polynomials_to_decompose), remove=True)

            # If there is a nontrivial primary decomposition, remove the integration variables from the `remainder_expression`
            for i,var in enumerate(all_integration_variables):
                if var not in integration_variables:
                    this_primary_sector_remainder_expression = this_primary_sector_remainder_expression.replace(i,1,remove=True)
                    break

            if use_symmetries and not split:
                # search for symmetries throughout the secondary decomposition
                indices = range(len(integration_variables))
                secondary_sectors = []
                for primary_sector in primary_sectors:
                    secondary_sectors.extend( strategy['secondary'](primary_sector, indices) )
                secondary_sectors = _reduce_sectors_by_symmetries\
                (
                    secondary_sectors,
                    'total number sectors',
                    indices,
                    use_iterative_sort,
                    use_light_Pak,
                    use_Pak,
                    dreadnaut_executable if use_dreadnaut else False,
                    name
                )
            else:
                secondary_sectors = strategy['secondary'](primary_sector, range(len(integration_variables)))

            # process the `secondary_sectors` in parallel
            if pool is not None:
                lowest_orders_and_function_declarations_and_pole_structures = \
                    pool.map(
                                _process_secondary_sector,
                                _make_environment( locals() )
                            )
            else:
                lowest_orders_and_function_declarations_and_pole_structures = \
                    list(map(
                            _process_secondary_sector,
                            _make_environment( locals() )
                        ))

            # get the `sector_index` after processing the secondary sectors
            sector_index = sector_index + len(lowest_orders_and_function_declarations_and_pole_structures)

            # update the global `lowest_orders`
            lowest_orders = np.min([item[0] for item in lowest_orders_and_function_declarations_and_pole_structures],axis=0)

            # update the global `function_declarations` and `pole_structures`
            for item in lowest_orders_and_function_declarations_and_pole_structures:
                _,f,p = item
                function_declarations.update(f)
                pole_structures.append(p)

    finally:
        # make sure the pool is closed
        if pool is not None:
            pool.close()

    # expand the `prefactor` to the required orders
    required_prefactor_orders = requested_orders - lowest_orders
    print('expanding the prefactor', prefactor, '(regulators:', regulators, ', orders:', required_prefactor_orders, ')')
    expanded_prefactor = expand_sympy(prefactor, regulators, required_prefactor_orders)

    # pack the `prefactor` into c++ function that returns a nested `Series`
    # and takes the `real_parameters` and the `complex_parameters`
    prefactor_type = 'secdecutil::Series<' * len(regulators) + 'integrand_return_t' + '>' * len(regulators)
    prefactor_function_body = _make_prefactor_function(expanded_prefactor, real_parameters, complex_parameters)

    # define the return type of "make_integrands"
    make_integrands_return_t = 'std::vector<' + 'secdecutil::Series<' * len(regulators) + \
                               'secdecutil::IntegrandContainer<integrand_return_t, real_t const * const' + \
                               '>' * (len(regulators) + 2)

    # generate translation from transform short names 'korobov#x#' and 'sidi#' to C++ macros
    pylink_qmc_instantiations_translation = generate_pylink_qmc_macro_dict('INSTANTIATE')
    pylink_qmc_extern_translation = generate_pylink_qmc_macro_dict('EXTERN')
    pylink_qmc_case_translation = generate_pylink_qmc_macro_dict('CASE')

    # parse the required pylink templates and generate a list of files to write
    pylink_qmc_transform_instantiation_rules = []
    pylink_qmc_extern_rules = []
    pylink_qmc_case_rules = []
    for pylink_qmc_transform in pylink_qmc_transforms:
        pylink_qmc_extern_rules.append(pylink_qmc_extern_translation[pylink_qmc_transform])
        pylink_qmc_case_rules.append(pylink_qmc_case_translation[pylink_qmc_transform])
    for i, pylink_qmc_transform in enumerate(chunks(pylink_qmc_transforms, 5)):
        pylink_qmc_transform_instantiation_rules.append(
            {
                'src': 'qmc_template_instantiations.cpp',
                'dest': 'qmc_template_instantiations_%i.cpp' % i,
                'replacements': {
                    'pylink_qmc_transforms': ' '.join([pylink_qmc_instantiations_translation[x] for x in pylink_qmc_transform])
                }
            }
        )

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
    template_replacements['integrand_getters'] = ''.join( 'nested_series_t<sector_container_t> get_integrand_of_sector_%i();\n' % i for i in range(1,sector_index+1) )
    template_replacements['sectors_initializer'] = ','.join( 'get_integrand_of_sector_%i()' % i for i in range(1,sector_index+1) )
    template_replacements['pole_structures_initializer'] = str(pole_structures).replace(' ','').replace("'","").replace('[','{').replace(']','}')
    template_replacements['pylink_qmc_externs'] = ' '.join(pylink_qmc_extern_rules)
    template_replacements['pylink_qmc_cases'] = ' '.join(pylink_qmc_case_rules)
    parse_template_file(os.path.join(template_sources, 'name.hpp'), # source
                        os.path.join(name,            name + '.hpp'), # dest
                        template_replacements)
    for filename in ['integrands.cpp', 'prefactor.cpp', 'pole_structures.cpp', 'functions.hpp']:
        parse_template_file(os.path.join(template_sources, 'src', filename),
                            os.path.join(name,             'src', filename),
                            template_replacements)
    for filename in ['pylink.cpp']:
        parse_template_file(os.path.join(template_sources, 'pylink', filename),
                            os.path.join(name,             'pylink', filename),
                            template_replacements)
    for pylink_qmc_transform in pylink_qmc_transform_instantiation_rules:
        #shared_keys = set(pylink_qmc_transform['replacements']).intersection(set(template_replacements))
        #if shared_keys:
        #    raise ValueError('`pylink_qmc_transforms` conflicts with `template_replacements`, rename "%s" in the `replacements` dictionary in `pylink_qmc_transforms`.' % ','.join(shared_keys) )
        pylink_qmc_transform_replacements = template_replacements.copy()
        pylink_qmc_transform_replacements.update(pylink_qmc_transform['replacements'])
        parse_template_file(os.path.join(template_sources, 'pylink', pylink_qmc_transform['src']),
                            os.path.join(name,             'pylink', pylink_qmc_transform['dest']),
                            pylink_qmc_transform_replacements)

    print('"' + name + '" done')

    # print message how to implement the dummy functions if applicable
    if have_dummy_functions:
        print(
                 "Declarations of the `functions` and their required derivatives are provided\n" + \
                 "in the file 'src/functions.hpp'. Please refer to that file for further\n" + \
                 "instructions."
             )

    # return the replacements in the template files
    return template_replacements


# namedtuple representing a make_package type package_generator (for use with sum_package)
MakePackage = namedtuple('MakePackage', list(inspect.signature(make_package).parameters))
# python <3.7 compatibility
MakePackage.__new__.__defaults__ = tuple([v.default for k,v in inspect.signature(make_package).parameters.items()
                                    if v.default != inspect._empty])
