#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1','k2'],
    external_momenta = ['p1','p2','p3','p4','p5'],
    #Lorentz_indices = ['mu'],

    propagators = ['k1**2','(k1-p1)**2','(k1-p1-p2)**2','(k1-p1-p2-p3)**2','k2**2','(k2-p1-p2-p3-p4)**2','(k1-k2)**2','(k1-k2+p4)**2','(k2-p1)**2','(k2-p1-p2)**2','(k2-p1-p2-p3)**2'],
    powerlist = [1,1,0,1,1,1,1,1,0,0,-1],

    #numerator = 'k1(mu)*k2(mu)',

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p1*p2', 'v1/2'),
                            ('p2*p3', 'v2/2'),
                            ('p1*p3', '(v4-v1-v2)/2'),
                            ('p1*p4', '(v2-v5-v4)/2'),
                            ('p2*p4', '(-v2-v3+v5)/2'),
                            ('p3*p4', 'v3/2')
                             ]
    )


    Mandelstam_symbols = ['v1', 'v2','v3', 'v4', 'v5']
    #mass_symbols = []


    loop_package(

    name = 'I73_29',

    loop_integral = li,

    real_parameters = Mandelstam_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    additional_prefactor = 'exp(-2*EulerGamma*eps)',

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '2G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',
    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = False,

    # no symmetries --> no need to run the full symmetry finder
    #use_Pak = False,

    )

    # generates 34 sectors, no symmetries
