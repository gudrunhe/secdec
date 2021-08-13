"""
Expansion by Regions
--------------------

Routines to perform an expansion by regions, see e.g. [PS11]_.

"""

from .polytope import Polytope
from .algebra import Polynomial, ExponentiatedPolynomial, Product
from .code_writer.make_package import MakePackage
from .misc import sympify_expression
import numpy as np

_sympy_one = sympify_expression(1)

def _poly_variable_transform(poly, new_variables, add_missing_var_to_coeff = True):
    """
    Transform the polynomial to use the new variables. The new variables can
    contain new variables and old variables in a different order, the
    polynomial will use them in the order given. If add_missing_var_to_coeff
    is True, the missing variables will be added to the coefficients, otherwise
    ignored (i.e. set to 1).
    """
    new_variables = list(map(sympify_expression,new_variables))

    extra_variables = [variable for variable in new_variables if variable not in poly.polysymbols]
    old_variables = poly.polysymbols + extra_variables
    missing_variables = [variable for variable in poly.polysymbols if variable not in new_variables]
    missing_variable_indexes = [old_variables.index(x) for x in missing_variables]
    newexpos = []
    newcoeffs = []

    order = [old_variables.index(x) for x in new_variables]

    for n, coeff in enumerate(poly.coeffs):
        if extra_variables:
            coeff_poly = Polynomial.from_expression(coeff, extra_variables)
        else:
            coeff_poly = Polynomial([[0]],[coeff])
        if add_missing_var_to_coeff:
            missing_variable_coeff = sympify_expression("1")
            for m in missing_variable_indexes:
                missing_variable_coeff *= old_variables[m]**poly.expolist[n][m]
            for k in range(len(coeff_poly.coeffs)):
                coeff_poly.coeffs[k] *= missing_variable_coeff
        newcoeffs.extend(coeff_poly.coeffs)
        for expo in coeff_poly.expolist:
            newexpos.append(np.append(poly.expolist[n],expo)[order])

    poly.coeffs = np.array(newcoeffs)
    poly.expolist = np.array(newexpos)
    poly.polysymbols = new_variables
    poly.number_of_variables = len(new_variables)

def _parse_expressions(polynomials_to_decompose, symbols_polynomials_to_decompose):
    """
    Return all the input polynomials converted to type ExponentiatedPolynomial with the
    variables symbols_polynomials_to_decompose. Also the exponent is a Polynomial with
    the same variables and all coefficients are sympy objects. Doesn't change the input
    polynomials.
    """
    parsed = []
    for n, poly in enumerate(list(polynomials_to_decompose)):
        if type(poly) is ExponentiatedPolynomial:
            poly = poly.copy()
            _poly_variable_transform(poly, symbols_polynomials_to_decompose)
            poly.exponent = Polynomial.from_expression(poly.exponent, symbols_polynomials_to_decompose)
            poly.coeffs = np.array([ sympify_expression(str(coeff)) for coeff in poly.coeffs ])
        elif type(poly) is Polynomial:
            poly = poly.copy()
            _poly_variable_transform(poly, symbols_polynomials_to_decompose)
            poly = ExponentiatedPolynomial(poly.expolist,
                                            np.array([ sympify_expression(str(coeff)) for coeff in poly.coeffs ]),
                                            Polynomial.from_expression(_sympy_one, symbols_polynomials_to_decompose), # exponent
                                            symbols_polynomials_to_decompose, copy=False)
        else:
            poly = sympify_expression(str(poly))
            if poly.is_Pow:
                assert len(poly.args) == 2
                poly_base = Polynomial.from_expression(poly.args[0], symbols_polynomials_to_decompose)
                poly_exponent = sympify_expression(str(poly.args[1]))
                poly = ExponentiatedPolynomial(poly_base.expolist,
                                                poly_base.coeffs,
                                                Polynomial.from_expression(poly_exponent, symbols_polynomials_to_decompose),
                                                symbols_polynomials_to_decompose, copy=False)
            else:
                poly = Polynomial.from_expression(poly, symbols_polynomials_to_decompose)
                poly = ExponentiatedPolynomial(poly.expolist,
                                                poly.coeffs,
                                                Polynomial.from_expression(_sympy_one, symbols_polynomials_to_decompose), # exponent
                                                symbols_polynomials_to_decompose, copy=False)
        parsed.append(poly)
    return parsed

def find_regions( exp_param_index , polynomial, indices = None, normaliz='normaliz', workdir='normaliz_tmp'):
    '''
    Find regions for the expansion by regions
    as described in [PS11]_.

    .. note::
        This function calls the command line executable of
        `normaliz` [BIR]_. See :ref:`installation_normaliz`
        for installation and a list of tested versions.

    :param exp_param_index:
        int;
        The index of the expansion parameter in the expolist.

    :param polynomials:
        an instance of :class:`.Polynomial`, for which to
        calculate the regions for.

    :param indices:
        list of integers or None;
        The indices of the parameters to be included
        in the asymptotic expansion. This should include all
        Feynman parameters (integration variables) and the
        expansion parameter. By default (``indices=None``),
        all parameters are considered.

    :param normaliz:
        string;
        The shell command to run `normaliz`.

    :param workdir:
        string;
        The directory for the communication with `normaliz`.
        A directory with the specified name will be created
        in the current working directory. If the specified
        directory name already exists, an :class:`OSError`
        is raised.

        .. note::
            The communication with `normaliz` is done via
            files.

    '''

    polytope_vertices = polynomial.expolist


    if indices is not None:
        dim = len(polytope_vertices[0])
        indices = list(indices)
        exp_param_index = indices.index(range(dim)[exp_param_index])
        polytope_vertices = [[vertex[i] for i in indices] for vertex in polytope_vertices]

    polytope = Polytope(vertices=polytope_vertices)
    polytope.complete_representation(normaliz, workdir)

    facets = polytope.facets[:,:-1] # do not need offset term "a_F"

    regions = facets[ facets[:,exp_param_index] > 0 ]

    if indices is not None:
        regions = np.array([[next(region) if i in indices else 0 for i in range(dim)] for region in [iter(r) for r in regions]])

    return regions

def apply_region(polynomials, region_vector, expansion_parameter_index):
    r'''
    Apply the `region_vector` to the input `polynomials`.

    .. note::
        Redefines the `expansion_parameter` as
        :math:`\rho \rightarrow \rho ^ n`, where
        :math:`n` is given by the `region_vector`.

    .. note::
        `apply_region` modifies the input
        `polynomials`.

    :param polynomials:
        iterable of polynomials;
        Polynomials to be computed in different regions.

    :param region_vector:
        vector-like array;
        Region specified by the power of the `expansion_parameter`
        in the rescaled variables. The region vectors have to be
        specified in the same order as the symbols are specified
        in the polynomials. E.g. if symbols are specified as
        ['x0','x1','rho'] and want rescaling x0 -> rho^i * x0,
        x1 -> rho^k * x1 and rho -> rho^n, then the region vector
        needs to be [i,k,n]

    :param expansion_parameter_index:
        integer;
        Index of the expansion parameter in the list of symbols.

    '''
    for polynomial in polynomials:
        polynomial.expolist[:,expansion_parameter_index] = np.dot(polynomial.expolist, region_vector)

    return polynomials


def derive_prod(poly_list,numerator,index,polynomial_name_indices):
    r"""
    Calculates the derivative of a product of polynomials using

    .. math::
        \frac{\partial}{\partial x_i} \left( \prod_j P_j^{\alpha_j} N \right)
        = \prod_j P_j^{\alpha_j-1} N'

    where N' is given by

    .. math::
        N' = \left( \sum_j N \alpha_j \frac{\partial P_j}{\partial x_i} \prod_{k \neq j} P_k \right)
        + \left(\prod_l P_l \right) \left[ \left( \sum_k \frac{\partial N}{\partial P_k}
        \frac{\partial P_k}{\partial x_i} \right) + \frac{\partial N}{\partial x_i} \right] .

    :param poly_list:
        list of :class:`.ExponentiatedPolynomial`;
        The exponentiated polynomials that should be differentiated.
        They need to be defined in terms of the symbols ``x0,x1,x2,..`` and
        ``p0,p1,p2..`` where ``p0,p1,p2..`` are the bases of the exponentiated
        polynomials.

    :param numerator:
        :class:`.Polynomial`;
        The numerator also defined as an exponentiated polynomial with
        ``symbols = [x0,x1,...,p0,p1,...]``.

    :param index:
        integer;
        Index of variable with respect to which the derivative is taken.

    :param polynomial_name_indices:
        iterable;
        Indices of polynomial names in poly_symbols.


    """
    number_of_polys = len(poly_list)
    symbols = numerator.polysymbols
    number_of_symbols = len(numerator.polysymbols)

    explicit_poly_list = []
    dummy_poly_list = []
    for i in range(number_of_polys):
        # define the polynomials without exponents because needed for explicit derivative
        explicit_poly = Polynomial(poly_list[i].expolist, poly_list[i].coeffs,poly_list[i].polysymbols)
        explicit_poly_list.append(explicit_poly)

        # define the polynomials as just 1*p1, 1*p2 etc. as only want to evaluate derivatives explicitly
        expolist = np.zeros(number_of_symbols,dtype=np.int)
        expolist[polynomial_name_indices[i]] = 1
        dummy_poly = Polynomial([expolist], [1], symbols)
        dummy_poly_list.append(dummy_poly)

    # calculate the product of all polynomials
    prod_all_polys = np.prod(dummy_poly_list)

    # calculate the new numerator according to formula:
    numerator_new = numerator.derive(index)*prod_all_polys

    for i in range(number_of_polys):
        prod = 1
        for j in range(number_of_polys):

            if j != i:

                prod *= dummy_poly_list[j]

        summand1 = numerator * poly_list[i].exponent * explicit_poly_list[i].derive(index) * prod
        summand2 = numerator.derive(polynomial_name_indices[i]) * explicit_poly_list[i].derive(index) * prod_all_polys

        numerator_new += (summand1+summand2)

    new_poly_list=[]
    for i in range(number_of_polys):
        new_poly_list.append(ExponentiatedPolynomial(poly_list[i].expolist, poly_list[i].coeffs,poly_list[i].exponent -1 ,poly_list[i].polysymbols))

    return new_poly_list,numerator_new

def expand_region(poly_list,numerator,index,order,polynomial_name_indices):
    r"""
    Expands the product of the polynomials in `poly_list` and the numerator with respect
    to the variable whose `index` is given to a desired order specified by `order`.

    :param poly_list:
        list of :class:`.ExponentiatedPolynomial`;
        The exponentiated polynomials that should be expanded.
        They need to be defined in terms of the symbols ``x0,x1,x2,..`` and
        ``p0,p1,p2..`` where ``p0,p1,p2..`` are the bases of the exponentiated
        polynomials.

    :param numerator:
        :class:`.Polynomial`;
        The numerator also defined as an exponentiated polynomial with
        ``symbols = [x0,x1,...,p0,p1,...]``.

    :param index:
        integer;
        Index of variable with respect to which the polynomials are expanded.

    :param order:
        integer;
        Desired order of expansion.

    :param polynomial_name_indices:
        list of int;
        Indices of polynomials in the symbols of the input polynomials.

    """

    if order < 0:
        return

    else:
        i = 1
        # 0th order term
        term = Product(*poly_list)
        term.factors.append(numerator)
        term = term.replace(index, 0, remove=True)
        for factor in term.factors:
            factor.simplify()
        yield term

        # calculate higher order terms
        derivative = derive_prod(poly_list,numerator,index,polynomial_name_indices)

        inverse_i_factorial = sympify_expression(1)
        while i <= order:

            inverse_i_factorial /= i
            term = Product(*derivative[0])
            term.factors.append((derivative[1] * inverse_i_factorial))
            term = term.replace(index, 0, remove=True)
            for factor in term.factors:
                if type(factor) is ExponentiatedPolynomial:
                    # make sure we don't simplify poly**0 to 1
                    if hasattr(factor.exponent, "simplify"):
                        factor.exponent = factor.exponent.simplify()
                    if type(factor.exponent) is not Polynomial and sympify_expression(factor.exponent).simplify() or \
                        type(factor.exponent) is Polynomial and (factor.exponent.expolist != 0).any():
                        factor.simplify()
                else:
                    factor.simplify()
            yield term

            derivative = derive_prod(derivative[0],derivative[1],index, polynomial_name_indices)

            i+=1


def make_regions(name, integration_variables, regulators, requested_orders, smallness_parameter,
                 polynomials_to_decompose, expansion_by_regions_order=0, real_parameters=[], complex_parameters=[],
                 normaliz='normaliz', polytope_from_sum_of = None, **make_package_args):
    r'''
    Applies the expansion by regions method
    (see e.g. [PS11]_) to a list of polynomials.

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
        The regulators of the integral.

    :param requested_orders:
        iterable of integers;
        Compute the expansion in the regulators to these
        orders.

    :param smallness_parameter:
        string or sympy symbol;
        The symbol of the variable in which the
        expression is expanded.

    :param polynomials_to_decompose:
        iterable of strings or sympy expressions or
        :class:`pySecDec.algebra.ExponentiatedPolynomial`
        or :class:`pySecDec.algebra.Polynomial`;
        The polynomials to be decomposed.

    :param expansion_by_regions_order:
        integer;
        The order up to which the expression is expanded
        in the `smallness_parameter`.
        Default: 0

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as real variables.

    :param complex_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as complex variables.

    :param normaliz:
        string;
        The shell command to run `normaliz`.
        Default: 'normaliz'

    :param polytope_from_sum_of:
        iterable of integers;
        If this value is None, the product of all the polynomials
        `polynomials_to_decompose` is used to determine the Newton
        polytope and the normal vectors to it. Otherwise, the sum of
        the polynomials with the indices given by this parameter are used.

    :param make_package_args:
        The arguments to be forwarded to
        :func:`pySecDec.code_writer.make_package`.

    '''
    # check that there are as many polynomial names as polynomials
    # and introduce more names if needed
    polynomial_names = []
    if 'polynomial_names' not in make_package_args or not make_package_args['polynomial_names']:
        polynomial_names = ['SecDecDecoPoly%i' % i for i in range(len(polynomials_to_decompose))]
    elif len(make_package_args['polynomial_names']) < len(polynomials_to_decompose):
        polynomial_names = make_package_args['polynomial_names'] + ['SecDecDecoPoly%i' % i for i in range(len(polynomials_to_decompose)-len( make_package_args['polynomial_names']))]
    else:
        polynomial_names = make_package_args['polynomial_names']

    symbols_polynomials_to_decompose = list(map(sympify_expression, integration_variables + regulators + polynomial_names + [smallness_parameter]))
    polynomial_name_indices = np.arange(len(polynomial_names)) + ( len(integration_variables) + len(regulators) )
    smallness_parameter_index = -1
    # only integration variables and the expansion parameter should be included in the calculation of the regions
    region_variable_indices = list(range(len(integration_variables))) + [ len(symbols_polynomials_to_decompose) - 1 ]

    # convert polynomials to ExponentedPolynomial with the proper variables
    polynomials_to_decompose = _parse_expressions(polynomials_to_decompose, symbols_polynomials_to_decompose)

    if not ( all( [np.count_nonzero(poly.expolist[:,len(integration_variables)+len(regulators):-1]) == 0 for poly in polynomials_to_decompose] ) ):
        raise NotImplementedError("Input polynomials for the asymptotic expansion are not allowed to depend on unexpanded polynomials.")

    # define the polynomial from which the Newton polytope is determined, for which we find the normal vectors
    if polytope_from_sum_of:
        polytope_polynomial = 0
        for index in polytope_from_sum_of:
            polytope_polynomial += Polynomial(polynomials_to_decompose[index].expolist, np.ones_like(polynomials_to_decompose[index].coeffs, dtype=int))
    else:
        polytope_polynomial = Polynomial(polynomials_to_decompose[0].expolist, np.ones_like(polynomials_to_decompose[0].coeffs, dtype=int))
        for poly in polynomials_to_decompose[1:]:
            polytope_polynomial *= Polynomial(poly.expolist, np.ones_like(poly.coeffs, dtype=int))

    # find the regions for the expansion (region vectors are always integer)
    regions = find_regions(smallness_parameter_index, polytope_polynomial, indices = region_variable_indices,
                            normaliz=normaliz, workdir=name)
    print("Total number of regions: {}".format(len(regions)))

    generators_args = []
    package_args_common = {'integration_variables':integration_variables, 'regulators':regulators,
                           'requested_orders':requested_orders, 'polynomial_names':polynomial_names}

    if regions.size > 0:
        for region in regions:
            # decompose the polynomials in the respective region
            polynomials_to_decompose_region_specific = apply_region([poly.copy() for poly in polynomials_to_decompose], region, smallness_parameter_index)

            # factor out the smallness_parameter and store its power
            polynomials_refactorized = []
            power_overall_smallness_parameter = 0
            for polynomial in polynomials_to_decompose_region_specific:
                factor0 ,factor1 = polynomial.refactorize(smallness_parameter_index).factors
                try:
                    exponent = factor0.exponent
                except AttributeError:
                    exponent = 1
                power_overall_smallness_parameter +=factor0.expolist[0][smallness_parameter_index] * exponent
                polynomials_refactorized.append(factor1)

            # compute overall power of smallness_parameter with all regulators -> 0
            power_overall_smallness_parameter_no_regulators = sympify_expression(power_overall_smallness_parameter)
            for regulator in regulators:
                power_overall_smallness_parameter_no_regulators = power_overall_smallness_parameter_no_regulators.subs(regulator,0)

            # define a dummy numerator
            numerator = Polynomial(np.zeros((1,len(symbols_polynomials_to_decompose)),dtype=int), np.array([1]), symbols_polynomials_to_decompose, copy=False)

            # exponent of the smallness parameter introduced by rescaling the integral measure
            power_smallness_parameter_measure = np.sum(region[:-1])

            print("Region ({}), smallness parameter lowest order is {}".format(
                ",".join(map(lambda x: str(sympify_expression(x)/region[smallness_parameter_index]),region[region_variable_indices])),
                (power_overall_smallness_parameter_no_regulators+power_smallness_parameter_measure)/region[smallness_parameter_index]))

            # expand the polynomials to the desired order:
            # expansion_by_regions_order - power_overall_smallness_parameter_no_regulators - power_smallness_parameter_measure
            expand_region_expansion_order = expansion_by_regions_order*region[smallness_parameter_index] - power_overall_smallness_parameter_no_regulators - power_smallness_parameter_measure
            series = expand_region(polynomials_refactorized, numerator, smallness_parameter_index, expand_region_expansion_order, polynomial_name_indices)

            for i,term in enumerate(series):
                package_args = make_package_args.copy()
                package_args['name'] = name + '_region_' + '_'.join(str(number).replace('-','m') for ind, number in enumerate(region) if ind in region_variable_indices) + '_expansion_order_' + str(power_overall_smallness_parameter_no_regulators + power_smallness_parameter_measure + i).replace('-','m').replace('/','d')

                package_args['other_polynomials'] = [term.factors.pop()]
                package_args['polynomials_to_decompose'] = term.factors
                prefactor = sympify_expression(smallness_parameter)**((power_overall_smallness_parameter+ power_smallness_parameter_measure +i)*(sympify_expression(1)/region[smallness_parameter_index]))
                if 'prefactor' in package_args:
                    prefactor *= sympify_expression(package_args['prefactor'])
                package_args['prefactor'] = prefactor

                package_args.update(package_args_common)

                generators_args.append(package_args)
    else:
        package_args = make_package_args.copy()
        package_args['name'] = name
        package_args['polynomials_to_decompose'] = map(str,polynomials_to_decompose)
        package_args.update(package_args_common)
        generators_args.append(package_args)

    return [MakePackage(**generator_args) for generator_args in generators_args]
