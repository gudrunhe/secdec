#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    banana = DistevalLibrary('banana_3mass/disteval/banana_3mass.json')

    # integrate
    str_result = banana(parameters={"s": 30.0, 'm1sq': 1.0, 'm2sq':1.3, 'm3sq':0.7}, verbose=True, epsrel = 1e-9, epsabs = 1e-100)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')