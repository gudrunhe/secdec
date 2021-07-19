#!/usr/bin/env python3

from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == '__main__':

    # example is case 8 of hep-ph/9704353, hep-ph/9605392

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['(k1+p1)**2','(k1-p2)**2','(k2+p1)**2','(k2-p2)**2','k2**2','(k1-k2)**2-msq'],
    loop_momenta = ['k1','k2'],
    external_momenta = ['p1','p2'],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p1*p2', 's/2')
                        ]
    )


    Mandelstam_symbols = ['s']
    mass_symbols = ['msq']

    loop_package(

    name = 'triangle2L_case8_full',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    additional_prefactor = '''gamma(1-2*eps)/(gamma(1+eps)*gamma(1+eps)*gamma(1-eps)*gamma(1-eps))*msq**(2+2*eps)''',

    # the highest order of the final epsilon expansion
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '1G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` and ``geometric``
    decomposition_method = 'geometric',
    # if you choose ``geometric`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    # whether or not to produce code to perform the contour deformation
    contour_deformation = True,

    )

    #li.exponentiated_U =(x4*x5 + x3*x5 + x2*x5 + x1*x5 + x1*x4 + x1*x3 + x1*x2 + x0*x5 + x0*x4 + x0*x3 + x0*x2)**(3*eps)

    #li.exponentiated_F =
    #( + (msq)*x4*x5**2 + (msq)*x3*x5**2 + (msq)*x2*x5**2 + (-s)*x2*x3*x5 + (msq)*x1*x5**2 + (msq)*x1*x4*x5 + (msq)*x1*x3*x5 + (msq - s)*x1*x2*x5 + (-s)*x1*x2*x3 + (msq)*x0*x5**2 + (msq)*x0*x4*x5 + (msq - s)*x0*x3*x5 + (msq)*x0*x2*x5 + (-s)*x0*x2*x3 + (-s)*x0*x1*x5 + (-s)*x0*x1*x4 + (-s)*x0*x1*x3 + (-s)*x0*x1*x2)**(-2*eps - 2)

