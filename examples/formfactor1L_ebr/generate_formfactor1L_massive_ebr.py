#!/usr/bin/env python3

from pySecDec.code_writer import sum_package, make_package
from pySecDec.make_regions import make_regions
import pySecDec as psd

if __name__ == "__main__":

    # Example 6 of Bernd Jantzen (arXiv:1111.2589)

    regions_generators = make_regions(

    # make_regions_args
    name = 'formfactor1L_massive_ebr',
    integration_variables = ['x0','x1','x2'],
    regulators = ['n','eps'],
    requested_orders = [0,0],
    smallness_parameter = 'z',
    polynomials_to_decompose = ['(x0+x1+x2)**((1+n/2)+(1+n/3)-3+2*eps)','((z*msq)*x2*(x0+x1+x2)+(qsq)*x0*x1)**(-(1+n/2)-(1+n/3)+1-eps)','x0**(n/2)','x1**(n/3)'],
    expansion_by_regions_order = 0,
    real_parameters = ['qsq','msq','z'],
    complex_parameters = [],

    # make_package_args
    polynomial_names = ['U','F'],
    contour_deformation_polynomial = 'F',
    positive_polynomials = ['U'],
    prefactor = '-gamma(eps + 1)', # limits n->0 already taken in prefactor
    decomposition_method = 'geometric',
    polytope_from_sum_of=[0,1])

    # sum_package
    sum_package('formfactor1L_massive_ebr',
        regions_generators, regulators = ['n','eps'], requested_orders = [0,0],
        real_parameters = ['qsq','msq','z'],
        complex_parameters = [])
