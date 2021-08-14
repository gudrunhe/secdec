#!/usr/bin/env python3
from pySecDec.loop_integral import LoopPackage
from pySecDec.code_writer.make_package import MakePackage
from pySecDec.code_writer.sum_package import Coefficient
from pySecDec.code_writer import sum_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['k1**2-msq','(k1+p1)**2-msq'],
    loop_momenta = ['k1'],
    replacement_rules = [
                            ('p1*p1', 's')
                        ]
    )
    Mandelstam_symbols = ['s']
    mass_symbols = ['msq']

    integral1 = LoopPackage(
        name = 'integral1',
        loop_integral = li,
        real_parameters = Mandelstam_symbols + mass_symbols,
        requested_orders = [0]
    )
    integral1_coeff = Coefficient(['s'],['msq'],['s','msq'])
    
    integral2 = MakePackage(
        name='integral2',
        integration_variables = ['z1','z2','z3'],
        regulators = ['eps'],
        polynomials_to_decompose = ['(s*z1+msq*z2*z3)**(-1+eps)'],
        polynomial_names = ['P'],
        contour_deformation_polynomial = 'P',
        real_parameters = Mandelstam_symbols + mass_symbols,
        requested_orders = [0]
    )
    integral2_coeff = Coefficient(['s**2'],['msq**2'],['s','msq'])

    sum_package(
        'mixed_packages',
        [integral1, integral2],
        regulators = ['eps'],
        coefficients = [[integral1_coeff, integral2_coeff]],
        real_parameters = Mandelstam_symbols + mass_symbols,
        requested_orders = [0]
    )
    
