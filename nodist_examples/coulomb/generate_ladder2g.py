#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

# integral with two gluon ladder rungs (box with 5p kinematics)
# 

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['msq',[1,2]],[0,[3,4]],[0,[6,7]],['msq',[2,6]],['msq',[1,7]],['msq',[3,6]],['msq',[4,7]]],
    #external_lines = [['p1',1],['p2',2],['p3',3],['p4+p5',4]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p45',4]],

    replacement_rules = [
                                ('p1*p1', 0),
                                ('p2*p2', 0),
                                ('p3*p3', 'msq'),
                                ('p45*p45', 's45'),
                                ('p1*p2', 's12/2'),
                                ('p1*p3', '(s45-s12-s23)/2'),
                                ('p1*p45', '(s23-s45)/2'),
                                ('p2*p45', '(-s12-s23+msq)/2'),
                                ('p2*p3', '(s23-msq)/2'),
                                ('p3*p45', '(s12-s45-msq)/2'),
                           ]
    )
    #replacement_rules = [       ('p45','p4+p5'),
    #                            ('p1*p1', 0),
    #                            ('p2*p2', 0),
    #                            ('p3*p3', 'msq'),
    #                            ('p4*p4', 'msq'),
    #                            ('p5*p5', 'mhsq'),
    #                            ('p1*p2', 's12/2'),
    #                            ('p1*p3', '(s45-s12-s23)/2'),
    #                            ('p1*p4', '(mhsq+s23-s51-s45)/2'),
    #                            ('p1*p5', '(s51-mhsq)/2'),
    #                            ('p2*p3', '(s23-msq)/2'),
    #                            ('p2*p4', '(-s23-s34+s51+msq)/2'),
    #                            ('p2*p5', '(s34-s12-s51)/2'),
    #                            ('p3*p4', 's34/2-msq'),
    #                            ('p3*p5', '(msq+s12-s34-s45)/2'),
    #                            ('p4*p5', '(s45-msq-mhsq)/2'),
    #                       ]


    Mandelstam_symbols = ['s12', 's23','s34','s45', 's51']
    mass_symbols = ['msq','mhsq']


    loop_package(

    name = 'ladder2g',

    loop_integral = li,

    real_parameters = Mandelstam_symbols+mass_symbols,
    #complex_parameters = mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,
    #requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 4,

    # the WorkSpace parameter for FORM
    form_work_space = '2G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    )
