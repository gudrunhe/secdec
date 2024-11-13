#! /usr/bin/env python
from pySecDec.code_writer import make_package
import pySecDec as psd
import sympy as sp

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromGraph(
        internal_lines = [[0,[3,4]],[0,[4,5]],[0,[3,6]],[0,[1,5]],[0,[1,6]],[0,[2,5]],[0,[2,6]]],
        external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
        replacement_rules = [
            ('p1*p1', 0),
            ('p2*p2', 0),
            ('p3*p3', 0),
            ('p4*p4', 'm2'),
            ('p1*p2', 's/2'),
            ('p2*p3', 't/2'),
            ('p1*p3', '(-s-t+m2)/2'),
            ('p3*p4', 's/2-m2/2'),
            ('p1*p4', 't/2-m2/2'),
            ('p2*p4', '(-s-t)/2')]
            )

    print(li.F)
    print(li.U)
    
    psd.loop_integral.draw.plot_diagram(
        internal_lines = [[0,[3,4]],[0,[4,5]],[0,[3,6]],[0,[1,5]],[0,[1,6]],[0,[2,5]],[0,[2,6]]],
        external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
        filename='BNP7_diagram')

    Mandelstam_symbols = ['s','t']
    mass_symbols = ['m2']

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

    print(Mandelstam_symbols+mass_symbols)

    make_package(
        name = 'BNP7',
        integration_variables = li.Feynman_parameters,
        regulators = ['eps'],
        requested_orders = [0],
        polynomials_to_decompose = [str(exponentiated_rescaled_U), str(exponentiated_rescaled_F)],
        polynomial_names = ['U','F'],
        other_polynomials = [jacobian_factor],
        prefactor = '('+additional_prefactor+')'+'*'+'('+str(li.Gamma_factor)+')',
        real_parameters = Mandelstam_symbols+mass_symbols,
        form_optimization_level = 2,
        form_work_space = '5G',
        contour_deformation_polynomial = 'F',
        positive_polynomials = 'U',
        decomposition_method = 'geometric')