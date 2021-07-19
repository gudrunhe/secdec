#!/usr/bin/env python3

from pySecDec.loop_integral import loop_package
import pySecDec as psd

# example is I3tilde of hep-ph/9605392

if __name__ == '__main__':

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['k1**2-msq','(k1+p1)**2-msq','(k1-k2)**2-msq','k2**2','(k2+p1)**2'],
    loop_momenta = ['k1','k2'],
    external_momenta = ['p1'],
    replacement_rules = [ ('p1*p1', 'psq') ]
    )


    Mandelstam_symbols = ['psq']
    mass_symbols = ['msq']

    loop_package(

    name = 'bubble2L_full',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    #additional_prefactor = '''gamma(1-2*eps)/(gamma(1+eps)*gamma(1+eps)*gamma(1-eps)*gamma(1-eps))*msq**(1+2*eps)''',
    # in hep-ph/9605392 the i^2 from i*Pi^(D/2) per loop is pulled out
    additional_prefactor = '-1',

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

    #li.exponentiated_U =(x2*x4 + x2*x3 + x1*x4 + x1*x3 + x1*x2 + x0*x4 + x0*x3 + x0*x2)**(3*eps - 1)

    #li.exponentiated_F =((-psq)*x2*x3*x4 + (msq)*x2**2*x4 + (msq)*x2**2*x3 + (-psq)*x1*x3*x4 + (2*msq)*x1*x2*x4 + (2*msq - psq)*x1*x2*x3 + (msq)*x1*x2**2 + (msq)*x1**2*x4 + (msq)*x1**2*x3 + (msq)*x1**2*x2 + (-psq)*x0*x3*x4 + (2*msq - psq)*x0*x2*x4 + (2*msq)*x0*x2*x3 + (msq)*x0*x2**2 + (2*msq - psq)*x0*x1*x4 + (2*msq - psq)*x0*x1*x3 + (2*msq - psq)*x0*x1*x2 + (msq)*x0**2*x4 + (msq)*x0**2*x3 + (msq)*x0**2*x2)**(-2*eps - 1)

