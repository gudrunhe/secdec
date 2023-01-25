#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    formfactor4L = DistevalLibrary('formfactor4L/disteval/formfactor4L.json')

    # integrate
    str_result = formfactor4L(parameters={"s": 90.0, 't': -2.5, 'pp4':1.6, 'msq':1.0}, verbose=True, epsrel = 1e-5, epsabs = 1e-100)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')