#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    Nbox2L = DistevalLibrary('Nbox2L_split_b/disteval/Nbox2L_split_b.json')

    # integrate
    str_result = Nbox2L(parameters={"s": -1.0, 't': -0.8, 'mm':0.1}, verbose=True, epsrel = 1e-4, epsabs = 1e-100)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')