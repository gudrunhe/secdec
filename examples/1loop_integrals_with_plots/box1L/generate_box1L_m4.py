#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
        internal_lines = [
            ['m1',[1,2]], 
            ['m2',[2,3]], 
            ['m3',[3,4]], 
            ['m4',[4,1]]
        ],
        external_lines = [
            ['p1',1], 
            ['p2',2], 
            ['p3',3], 
            ['p4',4],
        ],
        replacement_rules = [
            ('p4', '-p1-p2-p3'),
            ('p1*p1', 0),
            ('p2*p2', 0),
            ('p3*p3', 'p3sq'),
            ('p1*p2', 's/2'),
            ('p2*p3', 't/2-p3sq/2'),
            ('p1*p3', 'u/2-p3sq/2'),
            ('m1**2', 'm1sq'),
            ('m2**2', 'm2sq'),
            ('m3**2', 'm3sq'),
            ('m4**2', 'm4sq'),
        ]
    )

    psd.loop_package(
        name = 'box1L_m4',

        loop_integral = li,

        real_parameters = ['s','u','t','m1sq', 'm2sq', 'm3sq', 'm4sq', 'p3sq'],

        # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
        requested_orders = [0],
    )
