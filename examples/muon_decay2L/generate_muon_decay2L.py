#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
        internal_lines = [['mt',[1,4]],['mw',[4,2]],[0,[2,3]],[0,[4,5]], [0, [1,5]], ['mz', [5,3]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        replacement_rules = [
            ('p1*p1', 's'),
            ('p2*p2', 0),
            ('p3*p3', 0),
            ('p1*p2', 's/2'),
            ('p1*p3', 's/2'),
            ('p2*p3', 's/2'),
            ('mw**2', 'mwsq'),
            ('mz**2', 'mzsq'),
            ('mt**2', 'mtsq')
        ]
    )
    psd.loop_package(
        name = 'muon_decay2L',

        loop_integral = li,

        real_parameters = ['s', 'mwsq', 'mzsq', 'mtsq'],

        # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
        requested_orders = [0],

        # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
        form_optimization_level = 2,

        # the method to be used for the sector decomposition
        # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
        decomposition_method = 'geometric',
    )
