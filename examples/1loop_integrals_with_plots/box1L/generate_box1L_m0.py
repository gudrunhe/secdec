#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
        internal_lines = [
            [0, [1,2]], 
            [0, [2,3]], 
            [0, [3,4]], 
            [0, [4,1]]
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
        ]
    )

    psd.loop_package(
        name = 'box1L_m0',

        loop_integral = li,

        real_parameters = ['s','t','u','p3sq'],

        # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
        requested_orders = [0],
    )
