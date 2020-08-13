#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['m',[0,'a1']],['m',[0,'b1']],[0,['a1','b1']],['m',['a1','a2']],['m',['b1','b2']],[0,['a2','b2']],['m',['a2','a3']],['m',['b2','b3']],[0,['a3','b3']]],
    external_lines = [['p3',0],['p1','a3'],['p2','b3']],

    replacement_rules = [
                            ('p3*p3', 's'),
                            ('p1*p1', 'msq'),
                            ('p2*p2', 'msq'),
                            ('p1*p2', 's/2-msq'),
                            ('m**2', 'msq')
                        ]
    )


    Mandelstam_symbols = ['s']
    mass_symbols = ['msq']


    loop_package(

    name = 'ff3_massive',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '10G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric_ku',
    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    contour_deformation=True,

    split=True
    )
