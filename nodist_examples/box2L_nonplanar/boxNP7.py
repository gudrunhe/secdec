#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromGraph(

    internal_lines = [  [0,[1,5]],[0,[2,6]],[0,[1,2]],[0,[3,5]],[0,[3,6]],[0,[4,6]],[0,[4,5]]  ],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

    replacement_rules = [
                            ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '0'),
                            ('p4*p4', '0'),
                            ('p3*p2', 't/2'),
                            ('p1*p3', '-t/2-s/2'),
                            ('p1*p2', 's/2'),
                            ('p1*p4', 't/2'),
                            ('p2*p4', '-t/2-s/2'),
                            ('p3*p4', 's/2'),
                        ]

    )

    Mandelstam_symbols = ['s','t']
    mass_symbols = []


    loop_package(

    name = 'boxNP7',

    loop_integral = li,

    additional_prefactor='-gamma(3+2*eps)',

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '1G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` and ``geometric``
    decomposition_method = 'iterative',

    # there are singularities at one
    split = True,

    contour_deformation = True,

    )
