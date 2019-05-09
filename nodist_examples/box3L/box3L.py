#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    # Example used to demonstrate symmetry finder in Stephen Jones ACAT Proceedings 2017

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [ [0,[1,4]], [0,[1,5]], [0,[2,3]], [0,[2,7]], [0,[3,8]], [0,[4,6]], [0,[5,6]], [0,[5,7]], [0,[7,8]], [0,[6,8]] ],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
    powerlist=[1,1,1,1,1,1,1,1,1,1],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p1*p2', 's/2'),
                            ('p2*p3', 't/2'),
                            ('p1*p3', '-s/2-t/2')
                        ]

    )


    Mandelstam_symbols = ['s', 't']
    mass_symbols = []

    loop_package(
    name = 'box3L',
    loop_integral = li,
    real_parameters = Mandelstam_symbols + mass_symbols,
    requested_orders = [0],
    decomposition_method = 'iterative',
    contour_deformation = False,
    use_Pak = True,
    use_dreadnaut = False
    )

