#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(
        internal_lines = [['m',[1,2]],['m',[2,4]],['m',[4,6]],['m',[1,3]],['m',[3,5]],['m',[5,7]],[0,[2,3]],[0,[4,5]],[0,[6,7]]],
        external_lines = [['q',1],['p',6],['p',7]],
        replacement_rules = [
            ('q*q', 0),
            ('p*p', 'msq'),
            ('p*q', 0),
            ('m**2', 'msq')
        ]
    )

    Mandelstam_symbols = []
    mass_symbols = ['msq']

    psd.loop_package(
        name = 'g2_3L',

        loop_integral = li,

        real_parameters = Mandelstam_symbols + mass_symbols,

        # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
        requested_orders = [0],

        # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
        form_optimization_level = 2,

        # the method to be used for the sector decomposition
        # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
        decomposition_method = 'geometric',
        contour_deformation = False
    )
