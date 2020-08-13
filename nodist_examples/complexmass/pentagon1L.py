#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1'],
    external_momenta = ['p1','p2','p3','p4','p5'],
    #Lorentz_indices = ['mu'],

    propagators = ['k1**2','(k1+p1)**2-m1sq','(k1+p1+p2)**2-m2sq','(k1+p1+p2+p3)**2-m3sq','(k1+p1+p2+p3+p4)**2'],
    powerlist = [1,1,1,1,1],

    #numerator = 1,

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p5*p5', 0),
                       ('p1*p2', 's12/2'),
                       ('p1*p3', '(s45-s12-s23)/2'),
                       ('p1*p4', '(s23-s51-s45)/2'),
                       ('p1*p5', 's51/2'),
                       ('p2*p3', 's23/2'),
                       ('p2*p4', '(-s23-s34+s51)/2'),
                       ('p2*p5', '(s34-s12-s51)/2'),
                       ('p3*p4', 's34/2'),
                       ('p3*p5', '(s12-s34-s45)/2'),
                       ('p4*p5', 's45/2')
                                    ]

    )


    Mandelstam_symbols = ['s12','s23','s34','s45','s51']
    mass_symbols = ['m1sq','m2sq','m3sq']


    loop_package(

    name = 'pentagon_1L',

    loop_integral = li,

    real_parameters = Mandelstam_symbols,
    complex_parameters = mass_symbols,

    additional_prefactor='''gamma(1-2*eps)/(gamma(1+eps)*gamma(1-eps)**2)''',

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '500M',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    #decomposition_method = 'iterative',
    decomposition_method = 'geometric',
    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = True,

    # no symmetries --> no need to run the full symmetry finder
    #use_Pak = False,

    )

