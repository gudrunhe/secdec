import numpy as np
import sympy as sp
from pySecDec.algebra import Polynomial

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
                expansion_by_regions_order=0, normaliz='normaliz'):
    """
    Apply expansion of regions method to the loop integral using the function
    :func:`pySecDec.region_expand.make_regions`. See that function for the
    meaning of the parameters.

    :param loop_integral:
        :class:`pySecDec.loop_integral`;
        The loop integral to which the expansion by regions method is applied.

    """
    polynomials_to_decompose = [loop_integral.exponentiated_U, loop_integral.exponentiated_F] + loop_integral.measure.factors + [loop_integral.numerator]

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
        normaliz=normaliz,

        # make_package_args
        polynomial_names = polynomial_names,
        contour_deformation_polynomial = 'F',
        positive_polynomials = ['U'],
        prefactor = loop_integral.Gamma_factor,
        decomposition_method = 'iterative',
    )