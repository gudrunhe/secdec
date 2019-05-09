#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    # integral from 2-loop ttbar, called ggtt1 in 1303.1157

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['m1sq',[3,4]],[0,[4,5]],[0,[3,6]],['m2sq',[1,5]],['m2sq',[1,6]],['m2sq',[2,5]],['m2sq',[2,6]]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 'm1sq'),
                            ('p4*p4', 'm1sq'),
                            ('p1*p2', 's/2'),
                            ('p2*p3', 't/2-m1sq/2'),
                            ('p1*p3', '-s/2-t/2+m1sq/2'),
                            ('p3*p4', 's/2-m1sq'),
                            ('p1*p4', 't/2-m1sq/2'),
                            ('p2*p4', '-s/2-t/2+m1sq/2')
                        ]

    )


    Mandelstam_symbols = ['s','t']
    mass_symbols = ['m1sq','m2sq']


    loop_package(

    name = 'ggtt1',

    loop_integral = li,

    real_parameters = Mandelstam_symbols+mass_symbols,
    #complex_parameters = mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '1G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',
    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    )
