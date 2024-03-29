#!/usr/bin/env python3
import pySecDec as psd
from pySecDec.code_writer.make_package import make_package
from pySecDec.loop_integral import loop_package

if __name__ == "__main__":

    # 4-photon amplitude M++-- :
    #Amp4y_ppmm=-8*( 1 + (t**2+u**2)/s *Box6dim(t,u) +  (t-u)/s*(  BubbleD(u)-BubbleD(t) ) );

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [[0,[1,2]],[0,[2,1]]],
    external_lines = [['p1',1],['p2',2]],

    replacement_rules = [('p1*p1', 'uORt'),('p2*p2', 'uORt'),('p1*p2', 'uORt')]
    )

    Mandelstam_symbols = ['uORt']
    mass_symbols = []

    loop_package(

    name = 'yyyy_bubble',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` and ``geometric``
    decomposition_method = 'iterative',

    # pass code_writer.make_package to generate the old style pySecDec output
    package_generator=make_package

    )
