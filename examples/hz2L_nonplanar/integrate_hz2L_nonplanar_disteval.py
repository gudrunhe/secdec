#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    hz2L = DistevalLibrary('hz2L_nonplanar/disteval/hz2L_nonplanar.json')

    # integrate
    str_result = hz2L(parameters={"s12": 200.0, 's23': -23.0, 'mt2':9.0, 'mH2':1.56, 'mZ2':0.81}, verbose=True, epsrel = 1e-4, epsabs = 1e-100)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')