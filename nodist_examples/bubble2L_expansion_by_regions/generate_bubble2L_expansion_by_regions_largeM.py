from pySecDec.region_expand import make_regions
from pySecDec.algebra import Polynomial

import pySecDec as psd

regions_generator_args = make_regions(

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
decomposition_method = 'iterative'
)

psd.code_writer.sum_package('bubble2L_largeM',
    [psd.make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['eps'],requested_orders = [0],
    real_parameters = ['psq','msq','z'],
    complex_parameters = [],)

import configure