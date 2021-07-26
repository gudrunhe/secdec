#!/usr/bin/env python3

from pySecDec.code_writer.sum_package import sum_package, Coefficient
from pySecDec.loop_integral import LoopPackage, loop_package, LoopIntegralFromPropagators

if __name__ == "__main__":

    Mandelstam_symbols = ['s','t']
    mass_symbols = []

    real_parameters = Mandelstam_symbols + mass_symbols

    def make_1L_integral(powerlist, dim='4-2*eps', name=None):

        return LoopPackage(

            name = name if name else 'I_' + '_'.join( map(str,powerlist) ),

            loop_integral = LoopIntegralFromPropagators(

                propagators = ['(p1)**2','(p1+k1)**2','(p1+k1+k2)**2','(p1+k1+k2+k3)**2'],
                loop_momenta = ['p1'],
                external_momenta = ['k1','k2','k3','k4'],
                replacement_rules = [
                                    ('k1*k1', '0'),
                                    ('k2*k2', '0'),
                                    ('k3*k3', '0'),
                                    ('k1*k2', 's/2'),
                                    ('k1*k3', '-s/2-t/2'),
                                    ('k3*k2', 't/2'),
                                ],

                dimensionality = dim,

                powerlist = powerlist

            ),

            real_parameters = real_parameters,

            # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
            form_optimization_level = 2,

            # the WorkSpace parameter for FORM
            form_work_space = '100M',

            # the method to be used for the sector decomposition
            # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
            decomposition_method = 'iterative',
            # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
            # $PATH, you can set the path to the 'normaliz' command-line
            # executable here
            #normaliz_executable='/path/to/normaliz',

        )

    integrals = [
                    make_1L_integral([0,1,1,1]),
                    make_1L_integral([1,0,1,1]),
                    make_1L_integral([1,1,0,1]),
                    make_1L_integral([1,1,1,0]),
                    make_1L_integral([1,1,1,1], dim='6-2*eps'),
                ]

    coefficients = [
                       [
                            Coefficient(['1'], ['s'], ['eps'], real_parameters),
                            Coefficient(['1'], ['t'], ['eps'], real_parameters),
                            Coefficient(['1'], ['s'], ['eps'], real_parameters),
                            Coefficient(['1'], ['t'], ['eps'], real_parameters),
                            Coefficient(['2*(s+t)'], ['s*t'], ['eps'], real_parameters),
                       ]
                   ]

    # reduced
    sum_package(
        'D0_reduced',
        integrals,
        regulators = ['eps'],
        requested_orders = [0],
        coefficients = coefficients,
        real_parameters = real_parameters
    )

    # direct
    sum_package(
        'D0_direct',
        [make_1L_integral([1,1,1,1])],
        regulators = ['eps'],
        requested_orders = [0],
        real_parameters = real_parameters
    )

