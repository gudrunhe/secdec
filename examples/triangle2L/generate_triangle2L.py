#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
    internal_lines = [['m',[3,4]],['m',[4,5]],['m',[3,5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]],
    external_lines = [['p1',1],['p2',2],['p3',3]],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 's'),
                            ('p4*p4', 0),
                            ('p1*p2', 's/2'),
                            ('p2*p3', '-s/2'),
                            ('p1*p3', '-s/2'),
                            ('m**2', 'msq')
                        ]
    )


    Mandelstam_symbols = ['s']
    mass_symbols = ['msq']


    psd.loop_package(

    name = 'triangle2L',

    loop_integral = li,

    real_parameters = Mandelstam_symbols,
    complex_parameters = mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '100M',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    )
