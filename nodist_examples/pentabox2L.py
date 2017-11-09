#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd, sympy as sp

loop_momenta = k1, k2 = sp.symbols('k1 k2')
external_momenta = p1, p12, p123, p1234 = sp.symbols('p1 p12 p123 p1234')
Mandelstam_symbols = s12, s23, s34, s45, s51, x = sp.symbols('s12 s23 s34 s45 s51 x')
replacement_rules = [
                        (p1**2, 0), (p12**2, s12), (p123**2, s45), (p1234**2, 0),
                        (p1*p12, s12/2), (p1*p123, (s45 - s23)/2), (p1*p1234, -s51/2),
                        (p12*p123, (s12 + s45)/2), (p12*p1234, (s12 - s34)/2),
                        (p123*p1234, s45/2)
                    ]
propagators = [
                    k1**2, (k1 + x * p1)**2, (k1 + x * p12)**2, (k1 + p123)**2, (k1 + p1234)**2,
                    k2**2, (k2 - x * p1)**2, (k2 - x * p12)**2, (k2 - p123)**2, (k2 - p1234)**2,
                    (k1 + k2)**2
              ]

li = psd.loop_integral.LoopIntegralFromPropagators(
        propagators = propagators,
        loop_momenta = loop_momenta,
        external_momenta = external_momenta,
        replacement_rules = replacement_rules,
        powerlist = [1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1],
    )

loop_package(
    name = 'pentabox',
    loop_integral = li,
    use_Pak = False, # no symmetries
    real_parameters = Mandelstam_symbols,
    requested_order = 0,
    decomposition_method = 'geometric',
    form_work_space = '1G',
    # additional_prefactor = 'exp(2*EulerGamma*eps)',
)

# pentabox from arxiv:1511.09404
# analytic result available in ancillary files (Mathematica Notebook)

# s12 -> 3.1604938271604937, s23 -> -2.2115856637445197, s34 -> 0.9876543209876543, s45 -> 0.7901234567901234, s51 -> -1.8026785963650396, xp -> 1.125
# result (with additional prefator 'exp(2*EulerGamma*eps)'):
#    - 0.176721 * eps^-4
#    + (0.298436 -0.395075 I)*eps^-3
#    + (0.0865349 +1.13206 I)*eps^-2
#    + ((-2.6549-3.2401 I)+(0.0060391 -0.0121267 I) r1)*eps^-1
#    + (14.9462 -9.49785 I)+(0.0135741 -0.0114979*I)*r1
# with r1 = 0.912139*I
