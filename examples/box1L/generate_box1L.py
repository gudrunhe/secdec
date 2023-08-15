#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
        internal_lines = [['m',[1,2]], ['0',[2,3]], ['0',[3,4]], ['0',[4,1]]],
        external_lines = [['p1',1], ['p2',2], ['p3',3], ['p4',4]],
        replacement_rules = [
            ('p4', '-p1-p2-p3'),
            ('p1*p1', 's1'),
            ('p2*p2', 0),
            ('p3*p3', 0),
            ('p1*p2', 's/2-s1/2'),
            ('p1*p3', '-s/2-t/2'),
            ('p2*p3', 't/2'),
            ('m**2', 'msq')
        ]
    )

    Mandelstam_symbols = ['s','t','s1']
    mass_symbols = ['msq']

    psd.loop_package(
        name = 'box1L',

        loop_integral = li,

        real_parameters = Mandelstam_symbols + mass_symbols,

        # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
        requested_orders = [0],
    )
