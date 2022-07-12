#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(

    loop_momenta = ['k1','k2'],
    external_momenta = ['p1','p2','p3'],

    #INT(PL1x12,7,247,7,1,[],1,1,1,-1,1,1,1,1,0);
    propagators = ['k1**2','k2**2','(k1-k2)**2','(k1-p1)**2','(k2-p1)**2',
                   '(k1-p1-p2)**2','(k2-p1-p2)**2','(k1-p1-p2+p3)**2',
                   '(k2-p1-p2+p3)**2'],
    powerlist = [1,1,1,-1,1,1,1,1,0],

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 'm2'),
                            ('p1*p2', 's/2'),
                            ('p1*p3', '(m2-t)/2'),
                            ('p2*p3', '(s+t-m2)/2')
                        ]

    )

    Mandelstam_symbols = ['s', 't']
    mass_symbols = ['m2']


    loop_package(

    name = 'crossint5',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # additional_prefactor = '-((gamma(1-eps)**2*gamma(1+eps))/gamma(1-2*eps)/exp(-EulerGamma*eps)/3**eps)**(-2)/gamma(2*eps + 2)',

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [2],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '500M',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    )
