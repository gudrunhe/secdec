import numpy as np
import sympy as sp
from pySecDec.algebra import Polynomial, ExponentiatedPolynomial

def _poly_variable_transform(poly, new_variables, add_missing_var_to_coeff = True):
    """
    Transform the polynomial to use the new variables. The new variables can
    contain new variables and old variables in a different order, the
    polynomial will use them in the order given. If add_missing_var_to_coeff
    is True, the missing variables will be added to the coefficients, otherwise
    ingored (i.e. set to 1).
    """
    new_variables = list(map(sp.sympify,new_variables))

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
            missing_variable_coeff = sp.sympify("1")
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

def loop_regions(name, loop_integral, smallness_parameter,
                expansion_by_regions_order=0,
                contour_deformation=True,
                additional_prefactor=1, form_optimization_level=2,
                form_work_space='500M',
                add_monomial_regulator_power=None,
                decomposition_method='iterative',
                normaliz_executable='normaliz',
                enforce_complex=False,
                split=False, ibp_power_goal=-1,
                use_iterative_sort=True, use_light_Pak=True,
                use_dreadnaut=False, use_Pak=True,
                processes=None):
    """
    Apply expansion of regions method to the loop integral using the function
    :func:`pySecDec.region_expand.make_regions`. See that function for the
    meaning of the parameters.

    :param loop_integral:
        :class:`pySecDec.loop_integral`;
        The loop integral to which the expansion by regions method is applied.

    :param add_monomial_regulator_power:
        string or sympy symbol;
        Name of the regulator, using which monomial factors of the form
        $x_i**(n/p_i)$ are added, to regulate the integrals arising from
        the expansion by regions.

    """
    polynomials_to_decompose = [loop_integral.exponentiated_U, loop_integral.exponentiated_F] + loop_integral.measure.factors + [loop_integral.numerator]
    
    # add regulators of the form x_i**(n/p_i), where n is a regulator
    if add_monomial_regulator_power is not None:
        regulator = sp.sympify(add_monomial_regulator_power)
        loop_integral.regulators.insert(1,regulator)
        primes = [sp.prime(n+1) for n in range(len(loop_integral.integration_variables)-1)]
        monomial_factors = []
        for prime, variable in zip(primes,loop_integral.integration_variables[:-1]):
            variable = Polynomial.from_expression(variable,loop_integral.integration_variables)
            monomial = ExponentiatedPolynomial(variable.expolist,variable.coeffs,exponent=regulator/prime, polysymbols=variable.polysymbols)
            monomial_factors.append(monomial)
        variable = Polynomial.from_expression(loop_integral.integration_variables[-1],loop_integral.integration_variables)
        monomial_factors.append(ExponentiatedPolynomial(variable.expolist,variable.coeffs,exponent=sum(-regulator/prime for prime in primes), polysymbols=variable.polysymbols))
        polynomials_to_decompose = [loop_integral.exponentiated_U, loop_integral.exponentiated_F] + loop_integral.measure.factors + monomial_factors + [loop_integral.numerator]

    polynomial_names = ["U","F"] + ['Poly%i' % i for i in range(len(polynomials_to_decompose)-2)]

    variable_names = loop_integral.integration_variables + loop_integral.regulators + polynomial_names + [smallness_parameter]

    if "make_regions" not in dir():
        from ..region_expand import make_regions

    for poly in polynomials_to_decompose:
        _poly_variable_transform(poly, variable_names)

    return make_regions(
        # make_regions_args
        name = name,
        integration_variables = loop_integral.integration_variables,
        regulators = loop_integral.regulators,
        requested_orders = [],
        smallness_parameter = smallness_parameter,
        polynomials_to_decompose = polynomials_to_decompose,
        expansion_by_regions_order = expansion_by_regions_order,
        real_parameters = [],
        complex_parameters = [],
        polytope_from_sum_of = [0,1],
        normaliz=normaliz_executable,

        # make_package_args
        polynomial_names = polynomial_names,
        prefactor = loop_integral.Gamma_factor,

        contour_deformation_polynomial = 'F' if contour_deformation else None,
        positive_polynomials = ['U'],

        form_optimization_level = form_optimization_level,
        form_work_space = form_work_space,

        decomposition_method = decomposition_method,

        normaliz_executable = normaliz_executable,

        use_iterative_sort = use_iterative_sort,
        use_Pak = use_Pak,
        use_dreadnaut = use_dreadnaut,
        use_light_Pak = use_light_Pak,

        enforce_complex = enforce_complex,
        ibp_power_goal = ibp_power_goal,
        split = split,
        processes = processes
    )