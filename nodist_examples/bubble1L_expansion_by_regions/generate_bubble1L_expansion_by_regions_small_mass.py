from pySecDec.region_expand import make_regions
from pySecDec.algebra import Polynomial

import pySecDec as psd

make_regions(

# make_regions_args
name = 'bubble1L_expansion_by_regions_small_mass',
integration_variables = ['x1','x2'],
regulators = ['eps'],
requested_orders = [0],
smallness_parameter = 'z',
polynomials_to_decompose = ['(x1 + x2)**(2*eps-2)','(-s*x1*x2 + z*msq*(x1+x2)**2 )**(-eps)'],
expansion_by_regions_order = 5,
real_parameters = ['s','msq'],
complex_parameters = [],

# make_package_args
polynomial_names = ['U','F'],
contour_deformation_polynomial = 'F',
positive_polynomials = ['U'],
prefactor = 'gamma(2-1*(4-2*eps)/2)',
decomposition_method = 'iterative'
)

# Expected Result for s=1, msq=0.01
# + ((0,0) +/- (0,0))*eps^-2 + ((1,-1.38219e-09) +/- (1.93775e-07,2.76537e-07))*eps^-1 + ((1.53572,3.07812) +/- (7.28103e-07,9.42614e-07)) + O(eps)
