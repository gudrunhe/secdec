from pySecDec.region_expand import make_regions
import pySecDec as psd

# Example 6 of Bernd Jantzen (arXiv:1111.2589)
regions_generator_args = make_regions(

# make_regions_args
name = 'formfactor1L_massive_expansion_by_regions',
integration_variables = ['x0','x1','x2'],
regulators = ['n0','n1','eps'],
requested_orders = [0,0,0],
smallness_parameter = 'z',
polynomials_to_decompose = ['(x0+x1+x2)**((1+n0)+(1+n1)-3+2*eps)','((z**2*msq)*x2*(x0+x1+x2)+(qsq)*x0*x1)**(-(1+n0)-(1+n1)+1-eps)','x0**n0','x1**n1'],
expansion_by_regions_order = 0,
real_parameters = ['qsq','msq','z'],
complex_parameters = [],

# make_package_args
polynomial_names = ['U','F'],
contour_deformation_polynomial = 'F',
positive_polynomials = ['U'],
prefactor = '-gamma(eps + 1)', # limits n0,n1->0 already taken in the prefactor
decomposition_method = 'geometric'
)

psd.code_writer.sum_package('formfactor1L_massive_expansion_by_regions',
    [psd.make_package]*len(regions_generator_args),
    regions_generator_args, regulators = ['n0','n1','eps'],requested_orders = [0,0,0],
    real_parameters = ['qsq','msq','z'],
    complex_parameters = [],)

import configure

# Expected result from Eq(6.30) of arXiv:1111.2589
# N[-1/qsq*(1/2*Log[qsq/msq]^2 + Log[qsq/msq]*Log[1 - msq/qsq] - PolyLog[2, msq/qsq] + Pi^2/3) /. {qsq -> 100, msq -> 1}, 20]
# -0.13837355735881397785 + O(eps) + O(n1)
