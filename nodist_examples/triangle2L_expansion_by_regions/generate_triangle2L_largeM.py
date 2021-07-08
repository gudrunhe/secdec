#!/usr/bin/env python3

from pySecDec.code_writer import sum_package
from pySecDec.make_regions import make_regions
from pySecDec.algebra import Polynomial

if __name__ == '__main__':

    regions_generators = make_regions(

    # make_regions_args
    name = 'triangle2L_case8_largeM',
    integration_variables = ['x0','x1','x2','x3','x4','x5'],
    regulators = ['eps'],
    requested_orders = [0],
    smallness_parameter = 'z',
    polynomials_to_decompose = ['(x4*x5 + x3*x5 + x2*x5 + x1*x5 + x1*x4 + x1*x3 + x1*x2 + x0*x5 + x0*x4 + x0*x3 + x0*x2)**(3*eps)','( + (msq)*x4*x5**2 + (msq)*x3*x5**2 + (msq)*x2*x5**2 + (-z*s)*x2*x3*x5 + (msq)*x1*x5**2 + (msq)*x1*x4*x5 + (msq)*x1*x3*x5 + (msq - z*s)*x1*x2*x5 + (-z*s)*x1*x2*x3 + (msq)*x0*x5**2 + (msq)*x0*x4*x5 + (msq - z*s)*x0*x3*x5 + (msq)*x0*x2*x5 + (-z*s)*x0*x2*x3 + (-z*s)*x0*x1*x5 + (-z*s)*x0*x1*x4 + (-z*s)*x0*x1*x3 + (-z*s)*x0*x1*x2)**(-2*eps - 2)'],
    expansion_by_regions_order = -1,
    real_parameters = ['s','msq','z'],
    complex_parameters = [],

    # make_package_args
    form_work_space = '2G',
    polynomial_names = ['U','F'],
    contour_deformation_polynomial = 'F',
    positive_polynomials = ['U'],
    #prefactor = '''gamma(1-2*eps)/(gamma(1+eps)*gamma(1+eps)*gamma(1-eps)*gamma(1-eps))*msq**(2+2*eps)''',
    #prefactor = '''gamma(2+2*eps)''',
    decomposition_method = 'geometric',
    polytope_from_sum_of=[0,1]
    )

    sum_package('triangle2L_case8_largeM',
        regions_generators, regulators = ['eps'], requested_orders = [0],
        real_parameters = ['s','msq','z'],
        complex_parameters = [])
