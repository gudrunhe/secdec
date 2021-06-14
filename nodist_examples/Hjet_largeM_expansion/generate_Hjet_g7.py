from pySecDec.code_writer import sum_package, make_package
from pySecDec.make_regions import make_regions
from pySecDec.algebra import Polynomial

import pySecDec as psd

# analytic result from 1711.09875, eq. (II.68)

#def g7(s12,mh2,mt2,eps):
#      analytic = -(1/mt2)**(1+eps)*( 1/2 + (1+eps)/24*(s12+mh2)/mt2 )
#      return analytic

# in eps-expanded form to order eps:

#def g7_0(s12,mh2,mt2):
#      analytic0 = -(mh2 + 12*mt2 + s12)/(24*mt2**2)
#      return analytic0

#def g7_eps(s12,mh2,mt2):
#      analytic1 =  -(mh2 + s12 + (mh2 + 12*mt2 + s12)*log(1/mt2))/(24*mt2**2)
#      return analytic1

regions_generator_args = make_regions(

name = 'Hjet_g7',
integration_variables = ['x1','x2','x3'],
smallness_parameter = 'zz',
regulators = ['eps'],
requested_orders = [1],


polynomial_names = ['U','F'],
contour_deformation_polynomial = 'F',
positive_polynomials = ['U'],
#prefactor = '-gamma(1+eps)',
prefactor = '-1',
polynomials_to_decompose = ['(x1+x2+x3)**(2*eps-1)','(zz*(-s12*x2*x3-mh2*x1*x3) + bmt2*(x1+x2+x3)**2)**(-1-eps)'],
expansion_by_regions_order = 1,
decomposition_method = 'iterative',
real_parameters = ['s12','mh2','bmt2','zz'],
complex_parameters = [],
polytope_from_sum_of=[0,1]
)

sum_package('Hjet_g7',
    [make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['eps'],requested_orders = [1],
    real_parameters = ['s12','mh2','bmt2','zz'],
    complex_parameters = [])
