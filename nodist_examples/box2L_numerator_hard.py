#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1','k2'],
    external_momenta = ['p1','p2','p3','p4'],
    Lorentz_indices = ['mu1','mu2','mu3','mu4'],

    propagators = ['k1**2','k2**2-msq','(k1+p1)**2','(k1+p1+p2)**2','(k1+k2+p1+p2)**2-msq','(k1+k2+p1+p2+p3)**2-msq','(k2+p1+p2+p3)**2-msq','(k2+p1)**2','(k2+p1+p2)**2'],
    powerlist = [1,1,1,1,1,1,1,0,0],

    numerator = 'k1(mu1)*p2(mu1)*k1(mu2)*p2(mu2)*k2(mu3)*p2(mu3)*k2(mu4)*p2(mu4)',

    replacement_rules = [
                            ('p1*p1',0),
                            ('p2*p2',0),
                            ('p3*p3','mhsq'),
                            ('p4*p4','mhsq'),
                            ('p1*p2','s/2'),
                            ('p1*p3','-t/2-s/2+mhsq/2'),
                            ('p1*p4','t/2-mhsq/2'),
                            ('p2*p3','t/2-mhsq/2'),
                            ('p2*p4','-t/2-s/2+mhsq/2'),
                            ('p3*p4','s/2-mhsq')
                        ]

    )


    Mandelstam_symbols = ['s', 't']
    mass_symbols = ['mhsq','msq']


    loop_package(

    name = 'box2L_numerator_hard',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    additional_prefactor = 'exp(-2*euler*eps)',

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = -1,

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

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = True,

    )
