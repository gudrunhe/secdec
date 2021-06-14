from pySecDec.code_writer import sum_package, make_package
from pySecDec.make_regions import make_regions
import pySecDec as psd

# Expected Result for qsq=100, lsq=1, psq=1 (Compared with Long Chen)
# + ((0,0) +/- (2.93878e-19,0))*eps^-2 + ((-1.6865e-11,8.49456e-11) +/- (1.25602e-07,7.58591e-08))*eps^-1 + ((-0.244975,3.27442e-10) +/- (1.20291e-06,1.07642e-06)) + O(eps)

# Section 3.1 (Method of Regions for the Sudakov form factor) of Thomas Becher (arXiv:1803.04310)
regions_generator_args = make_regions(

# make_regions_args
name = 'formfactor1L_massless_expansion_by_regions',
integration_variables = ['x0','x1','x2'],
regulators = ['eps'],
requested_orders = [0],
smallness_parameter = 'z',
polynomials_to_decompose = ['( + (1)*x0 + (1)*x1 + (1)*x2)**(2*eps - 1)','( + (z**2*lsq)*x1*x2 + (z**2*psq)*x0*x2 + (qsq)*x0*x1)**(-eps - 1)'],
expansion_by_regions_order = 0,
real_parameters = ['qsq','lsq','psq','z'],
complex_parameters = [],

# make_package_args
polynomial_names = ['U','F'],
contour_deformation_polynomial = 'F',
positive_polynomials = ['U'],
prefactor = '-gamma(eps + 1)', # non-singular limits 1+n0->1 and 1+n1->1 already taken in prefactor
decomposition_method = 'iterative',
polytope_from_sum_of=[0,1]
)

sum_package('formfactor1L_massless_expansion_by_regions',
    [make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['eps'],requested_orders = [0],
    real_parameters = ['qsq','lsq','psq','z'],
    complex_parameters = [],)
