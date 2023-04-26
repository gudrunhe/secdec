#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp
import numpy as np

if __name__ == "__main__":

    # load the library
    amplitude = DistevalLibrary('full_amplitude/disteval/full_amplitude.json')

    # integrate
    str_result = amplitude(parameters={"s": 3.0, "t": -1.0, "u": -2.0}, verbose=True, epsrel = 1e-10)
    
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    valueN =result[1].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")
    errorN = result[1].coeff("plusminus")

    #Express results in terms of fine structure constant (and add back common pi factors of coefficients)
    value *= np.pi**5 * 4**3
    error *= np.pi**5 * 4**3 
    valueN *= np.pi**5 * 4**3
    errorN *= np.pi**5 * 4**3

    # examples how to access individual orders
    print('Numerical Result Proportional to alpha**3 (Nf is the number of leptons):')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')', '\n', '        + Nf*(', valueN.coeff('eps',-2), '+/- (', errorN.coeff('eps',-2), '))', '\n')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')', '\n', '        + Nf*(', valueN.coeff('eps',-1), '+/- (', errorN.coeff('eps',-1), '))', '\n')
    print('eps^0:', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')', '\n', '       + Nf*(',valueN.coeff('eps',0), '+/- (', errorN.coeff('eps',0), '))')