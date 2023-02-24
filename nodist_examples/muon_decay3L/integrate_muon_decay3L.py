#!/usr/bin/env python3

import sympy as sp

if __name__ == "__main__":

    # load c++ library
    from pySecDec.integral_interface import DistevalLibrary
    name = 'muon_decay3L'
    loop_integral = DistevalLibrary('{0}/disteval/{0}.json'.format(name))

    str_result = loop_integral(parameters={'mwsq' : 0.78, 'mzsq' : 1.0}, verbose=True)