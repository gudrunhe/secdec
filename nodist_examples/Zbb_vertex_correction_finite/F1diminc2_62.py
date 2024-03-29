#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['k1^2', '(k1-k2)^2', '(k1-k2+p1)^2-mZ^2', '(k2)^2', '(k2+p2)^2', '(k1+p1+p2)^2', '(k2+p1)^2'],
    loop_momenta = ['k1','k2'],

    external_momenta = ['p1','p2','p3','p4'],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 's'),
                            ('p1*p2', 's/2'),
                            ('p2*p3', '-s/2'),
                            ('p1*p3', '-s/2'),
                            ('s', 'mZ^2'),
                            ('mZ', 1)
                        ],
    dimensionality='6-2*eps',
    powerlist=[0, 1, 2, 1, 1, 2, 0]
    )


    Mandelstam_symbols = []
    mass_symbols = []


    loop_package(

    name = 'F1diminc2_62',

    additional_prefactor = '-exp(2*EulerGamma*eps)',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [1],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '1G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'iterative',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = True
    )
