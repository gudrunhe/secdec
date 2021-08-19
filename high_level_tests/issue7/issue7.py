#! /usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromGraph(
        internal_lines = [[0,[1,2]], [0,[1,2]], [0,[1,2]]],
        external_lines = [['p1',1],['p2',2]],

        powerlist = [ 1, 1, 1],

        replacement_rules = [
            ('p1*p1', 's'),
            ('p2*p2', 's')
        ]
        ,
        dimensionality= '3-2*eps'
    )

    loop_package(
        
        name = 'issue7',
        loop_integral = li,
        real_parameters = ['s'],
        requested_order = -1,
        form_optimization_level = 2,
        form_work_space = '100M',
        contour_deformation = False,
        decomposition_method = 'iterative'
    )
