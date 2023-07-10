#! /usr/bin/env python
#from pySecDec.loop_integral import loop_package
from pySecDec.code_writer import make_package
import pySecDec as psd
import sympy as sp

if __name__ == "__main__":

    # massless non-planar 6-propagator box
    name = 'BNP6_wu'
    internal_lines = [[0,[1,2]],[0,[2,5]],[0,[3,5]],[0,[3,1]],[0,[4,5]],[0,[4,1]]]
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]]

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = internal_lines,
    external_lines = external_lines,

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p1*p2', 's/2'),
                            ('p2*p3', 't/2'),
                            ('p1*p3', '(-s-t)/2'),
                            ('p3*p4', 's/2'),
                            ('p1*p4', 't/2'),
                            ('p2*p4', '(-s-t)/2')
                        ]
    )
    #psd.loop_integral.draw.plot_diagram(internal_lines, external_lines,name)

    Mandelstam_symbols = ['s','t']
    mass_symbols = []

    # rescale xi -> ai * xi
    rescaled_U = sp.sympify(str(li.U))
    rescaled_F = sp.sympify(str(li.F))
    rescale_factors = []
    for i, fparam in enumerate(li.Feynman_parameters):
        rescaled_U = rescaled_U.subs(fparam,f'(a{i}*{fparam})')
        rescaled_F = rescaled_F.subs(fparam,f'(a{i}*{fparam})')
        Mandelstam_symbols.append(f'a{i}')
        rescale_factors.append(f'a{i}')
    jacobian_factor = '*'.join(rescale_factors)

    rescaled_U = psd.algebra.Polynomial.from_expression(str(rescaled_U),li.exponentiated_U.polysymbols)
    rescaled_F = psd.algebra.Polynomial.from_expression(str(rescaled_F),li.exponentiated_F.polysymbols)
    exponentiated_rescaled_U = psd.algebra.ExponentiatedPolynomial(rescaled_U.expolist,rescaled_U.coeffs,li.exponent_U,li.exponentiated_U.polysymbols)
    exponentiated_rescaled_F = psd.algebra.ExponentiatedPolynomial(rescaled_F.expolist,rescaled_F.coeffs,li.exponent_F,li.exponentiated_F.polysymbols)

    #additional_prefactor = '3*Gamma[1+2*eps]*Gamma[1-eps]^3/Gamma[1-3*eps]/(1+4*eps)'
    additional_prefactor = '1'

    make_package(

    name = name,

    integration_variables = li.Feynman_parameters,
    regulators = ['eps'],

    requested_orders = [0],

    polynomials_to_decompose = [str(exponentiated_rescaled_U), str(exponentiated_rescaled_F)],
    polynomial_names = ['U','F'],
    other_polynomials = [jacobian_factor],

    prefactor = '('+additional_prefactor+')'+'*'+'('+str(li.Gamma_factor)+')',
    real_parameters = Mandelstam_symbols+mass_symbols,

    form_optimization_level = 4,
    form_work_space = '5G',

    contour_deformation_polynomial = 'F',
    positive_polynomials = 'U',
    decomposition_method = 'geometric',

    )

