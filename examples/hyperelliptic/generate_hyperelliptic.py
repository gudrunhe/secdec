#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromPropagators(

    loop_momenta = ['k1','k2'],
    external_momenta = ['p1','p2','p3','p4'],

    propagators = ['k1**2-m1sq','(k1-p1)**2-m1sq','(k1-p1-p2)**2-m1sq','k2**2-m3sq','(k2-p3)**2-m3sq','(k1-k2-p1-p2-p3)**2-m2sq','(k1+k2)**2-m2sq'],
    powerlist = [1,1,1,1,1,1,1],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p1*p2', 's/2'),
                            ('p2*p3', 't/2'),
                            ('p1*p3', '-s/2') ]

    )


    Mandelstam_symbols = ['s', 't']
    mass_symbols = ['m1sq','m2sq','m3sq']


    psd.loop_package(

    name = 'hyperelliptic',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 4,

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = True,

    # no symmetries --> no need to run the full symmetry finder
    use_Pak = False,

    )

    # generates 34 sectors, no symmetries
