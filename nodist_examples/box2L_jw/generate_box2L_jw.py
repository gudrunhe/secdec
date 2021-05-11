#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['q2**2','q1**2-msq','(q2+k2)**2','(q1-k1-k2+k4)**2','(q1+q2-k1+k4)**2','(q1-k1)**2-msq','(q1+q2-k1)**2-msq','(q1+k2)**2','(q1+q2)**2'],
    powerlist = [1,1,1,1,1,2,1,0,0],
    loop_momenta = ['q1','q2'],
    replacement_rules = [
                            ('k1*k1',0),
                            ('k2*k2',0),
                            ('k4*k4','msq'),
                            ('k1*k2','s/2'),
                            ('k1*k4','-3*msq/28+s/2+t/2'),
                            ('k2*k4','msq/2-t/2')
                        ]
    )

    Mandelstam_symbols = ['s','t']
    mass_symbols = ['msq']


    loop_package(

    name = 'box2L_jw',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    #additional_prefactor = '(s/msq)**(3/2)',

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '2000M',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` and ``geometric``
    decomposition_method = 'iterative',
    # if you choose ``geometric`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    contour_deformation = True

    )
