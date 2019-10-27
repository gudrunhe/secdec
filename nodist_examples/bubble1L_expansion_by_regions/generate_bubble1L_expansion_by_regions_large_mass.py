from pySecDec.region_expand import make_regions
from pySecDec.algebra import Polynomial

import pySecDec as psd

regions_generator_args = make_regions(

# make_regions_args
name = 'bubble1L_expansion_by_regions_large_mass',
integration_variables = ['x1','x2'],
regulators = ['eps'],
requested_orders = [0],
smallness_parameter = 'z',
polynomials_to_decompose = ['(x1 + x2)**(2*eps-2)','(-z*s*x1*x2 + msq*(x1+x2)**2 )**(-eps)'],
expansion_by_regions_order = 2,
real_parameters = ['s','msq','z'],
complex_parameters = [],

# make_package_args
polynomial_names = ['U','F'],
contour_deformation_polynomial = 'F',
positive_polynomials = ['U'],
prefactor = 'gamma(2-1*(4-2*eps)/2)',
decomposition_method = 'iterative'
)

psd.code_writer.sum_package('bubble1L_expansion_by_regions_large_mass',
    [psd.make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['eps'],requested_orders = [0],
    real_parameters = ['s','msq','z'],
    complex_parameters = [],)