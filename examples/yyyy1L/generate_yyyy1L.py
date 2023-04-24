#! /usr/bin/env python3
import pySecDec as psd

# 4-photon amplitude M++--

if __name__ == '__main__':

    integrals = [
        # one loop bubble (u)
        psd.LoopPackage(
            name = "bubble_u",
            loop_integral = psd.loop_integral.LoopIntegralFromGraph(
                internal_lines = [[0,[1,2]],[0,[2,1]]],
                external_lines = [['p1',1],['p2',2]],
                replacement_rules = [('p1*p1', 'u'),('p2*p2', 'u'),('p1*p2', 'u')])),
        # one loop bubble (t)
        psd.LoopPackage(
            name = 'bubble_t',
            loop_integral = psd.loop_integral.LoopIntegralFromGraph(
                internal_lines = [[0,[1,2]],[0,[2,1]]],
                external_lines = [['p1',1],['p2',2]],
                replacement_rules = [('p1*p1', 't'),('p2*p2', 't'),('p1*p2', 't')])),
        # one loop box (in 6 dimensions)
        psd.LoopPackage(
            name = "box_6",
            loop_integral = psd.loop_integral.LoopIntegralFromGraph(
                internal_lines = [['0',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
                external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
                dimensionality = '6-2*eps',
                replacement_rules = [
                    ('p1*p1', 0),
                    ('p2*p2', 0),
                    ('p3*p3', 0),
                    ('p4*p4', 0),
                    ('p3*p2', 'u/2'),
                    ('p1*p2', 't/2'),
                    ('p1*p4', 'u/2'),
                    ('p1*p3', '-u/2-t/2'),
                    ('p2*p4', '-u/2-t/2'),
                    ('p3*p4', 't/2')
                ])),
        # one loop box (in 8 dimensions)
        psd.LoopPackage(
            name = "box_8",
            loop_integral = psd.loop_integral.LoopIntegralFromGraph(
                internal_lines = [['0',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
                external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
                dimensionality= '8-2*eps',
                replacement_rules = [
                    ('p1*p1', 0),
                    ('p2*p2', 0),
                    ('p3*p3', 0),
                    ('p4*p4', 0),
                    ('p3*p2', 'u/2'),
                    ('p1*p2', 't/2'),
                    ('p1*p4', 'u/2'),
                    ('p1*p3', '-u/2-t/2'),
                    ('p2*p4', '-u/2-t/2'),
                    ('p3*p4', 't/2')
                ]))
    ]

    coefficients = {
        "M++--": [
            # bubble (u) coefficient
            '-8*(t-u)/(-u-t)',

            # bubble (t) coefficient
            '-8*(u-t)/(-u-t)',

            # box6 coefficient
            '-8*(t^2+u^2)/(-u-t)',

            # box8 coefficient
            '-8*3*(2*eps)'
        ]
    }

    # generate code sum of (int * coeff)
    psd.sum_package(
        'yyyy1L',
        integrals,
        coefficients = coefficients,
        regulators = ['eps'],
        requested_orders = [0],
        real_parameters = ['t', 'u'],
        complex_parameters = []
    )
