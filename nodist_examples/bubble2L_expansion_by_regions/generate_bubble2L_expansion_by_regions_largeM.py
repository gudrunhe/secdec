#!/usr/bin/env python3

from pySecDec.code_writer import make_package, sum_package
from pySecDec.make_regions import make_regions
from pySecDec.algebra import Polynomial

import pySecDec as psd

# example is I3tilde of hep-ph/9605392

if __name__ == '__main__':

    regions_generators = make_regions(

    # make_regions_args
    name = 'bubble2L_largeM',
    integration_variables = ['x0','x1','x2','x3','x4'],
    regulators = ['eps'],
    requested_orders = [0],
    smallness_parameter = 'z',
    polynomials_to_decompose = ['(x2*x4 + x2*x3 + x1*x4 + x1*x3 + x1*x2 + x0*x4 + x0*x3 + x0*x2)**(3*eps - 1)','((-z*psq)*x2*x3*x4 + (msq)*x2**2*x4 + (msq)*x2**2*x3 + (-z*psq)*x1*x3*x4 + (2*msq)*x1*x2*x4 + (2*msq - z*psq)*x1*x2*x3 + (msq)*x1*x2**2 + (msq)*x1**2*x4 + (msq)*x1**2*x3 + (msq)*x1**2*x2 + (-z*psq)*x0*x3*x4 + (2*msq - z*psq)*x0*x2*x4 + (2*msq)*x0*x2*x3 + (msq)*x0*x2**2 + (2*msq - z*psq)*x0*x1*x4 + (2*msq - z*psq)*x0*x1*x3 + (2*msq - z*psq)*x0*x1*x2 + (msq)*x0**2*x4 + (msq)*x0**2*x3 + (msq)*x0**2*x2)**(-2*eps - 1)'],
    expansion_by_regions_order = 1,
    real_parameters = ['psq','msq','z'],
    complex_parameters = [],

    # make_package_args
    polynomial_names = ['U','F'],
    contour_deformation_polynomial = 'F',
    positive_polynomials = ['U'],
    #prefactor = '-1',
    decomposition_method = 'iterative',
    polytope_from_sum_of=[0,1]
    )

    # generate code that will calculate the sum of all regions and all orders in
    # the smallness parameter
    sum_package('bubble2L_largeM',
        regions_generators, regulators = ['eps'],requested_orders = [0],
        real_parameters = ['psq','msq','z'],
        complex_parameters = [])
