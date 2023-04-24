#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromPropagators(

    loop_momenta = ['k1','k2'],
    external_momenta = ['p1','p2','p3','p4'],

    propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2','(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2','(k1+p3)**2'],
    powerlist = [1,1,1,1,1,1,1,-1],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p1*p2', 's/2'),
                            ('p2*p3', 't/2'),
                            ('p1*p3', '-s/2-t/2')
                        ]

    )


    Mandelstam_symbols = ['s', 't']
    mass_symbols = []


    psd.loop_package(

    name = 'box2L_invprop',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = False,

    # no symmetries --> no need to run the full symmetry finder
    use_Pak = False,

    )
