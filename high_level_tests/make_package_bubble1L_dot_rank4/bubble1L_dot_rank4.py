#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['k^2-m1^2','(k-p1)^2-m2^2'],
    loop_momenta = ['k'],
    external_momenta = ['p1','p2'],

    Lorentz_indices = ['mu','nu'],
    numerator = 'k(mu)*k(mu)*k(nu)*k(nu)',

    powerlist = [3,1],
    replacement_rules = [
                            ('p1*p1', 's1'),
                            ('p2*p2', 's2'),
                            ('(p1+p2)*(p1+p2)','sij')
                        ]
    )

    Mandelstam_symbols = ['s1','s2','sij']
    mass_symbols = ['m1','m2']

    loop_package(

    name = 'bubble1L_dot_rank4',

    loop_integral = li,


    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    complex_parameters =['s1','s2','sij'],
    real_parameters =['m1','m2'],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 3,

    # the WorkSpace parameter for FORM
    form_work_space = '500M',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = True,

    pylink_qmc_transforms=['baker','korobov3'],

    package_generator=psd.code_writer.make_package

    )
