#!/usr/bin/env python3

import sympy as sp

if __name__ == "__main__":

    # load c++ library
    from pySecDec.integral_interface import DistevalLibrary
    name = 'muon_decay2L'
    loop_integral = DistevalLibrary('{0}/disteval/{0}.json'.format(name))

    #integrate
    str_result = loop_integral(parameters={'s' : 3, 'mwsq' : 0.78, 'mzsq' : 1.0, 'mtsq' : 0.00038}, verbose=True)
