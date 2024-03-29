#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
    internal_lines = [[0,[1,2]],[0,[2,5]],[0,[1,6]],['mt',[3,5]],['mt',[3,6]],['mt',[4,5]],['mt',[4,6]]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 'mH2'),
                            ('p4*p4', 'mZ2'),
                            ('mt*mt', 'mt2'),
                            ('p1*p2', 's12/2'),
                            ('p2*p3', '(s23-mH2)/2'),
                            ('p1*p3', '(-s12-s23+mZ2)/2'),
                            ('p3*p4', '(s12-mH2-mZ2)/2'),
                            ('p1*p4', '(s23-mZ2)/2'),
                            ('p2*p4', '(-s12-s23+mH2)/2')
                        ]

    )


    Mandelstam_symbols = ['s12','s23']
    mass_symbols = ['mt2','mH2','mZ2']


    psd.loop_package(

    name = 'hz2L_nonplanar',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 4,

    # the WorkSpace parameter for FORM
    form_work_space = '200M',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    # there are no symmetries
    use_dreadnaut = False,
    use_Pak = False,

    )
