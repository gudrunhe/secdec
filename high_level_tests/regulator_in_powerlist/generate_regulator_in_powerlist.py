#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    #li = psd.loop_integral.LoopIntegralFromPropagators(
    #propagators = ['(k1-p1)**2','(k2-p1)**2','k2**2','k1**2','(k2-k1)**2'],
    #loop_momenta = ['k1','k2'],

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [[0,[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]],[0,[2,4]]],
    external_lines = [['p1',1],['p1',3]],

    powerlist = ['eps',1,'eps',1,1],
    dimensionality = '2-2*eps',
    regulator = 'eps',
    replacement_rules = [('p1*p1','s')]
    )

    loop_package(
    name = 'regulator_in_powerlist',
    loop_integral = li,
    requested_orders = [0],
    real_parameters = ['s'],
    decomposition_method = 'geometric',
    )
