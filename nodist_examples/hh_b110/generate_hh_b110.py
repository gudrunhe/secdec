#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(

    propagators = ['(l1)**2-mt2','(l1+q1p)**2-mt2','2*l1*q2m+s12'],
    powerlist = ['1','1','b2'],
    loop_momenta = ['l1'],
    replacement_rules = [
                            ('q1p*q1p','0'),
                            ('q2m*q2m','0'),
                            ('q3m*q3m','0'),
                            ('q1p*q2m','s12/2'),
                            ('q1p*q3m','-s13/2'),
                            ('q2m*q3m','0')
                        ],
    regulators = ['b2','eps']
    )

    Mandelstam_symbols = ['s12','s13']
    mass_symbols = ['mt2']

    loop_package(
    name = 'hh_b110',
    loop_integral = li,
    #additional_prefactor='-gamma(3+2*eps)',
    requested_orders = [0,0],
    real_parameters = Mandelstam_symbols + mass_symbols,
    )
