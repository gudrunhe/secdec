#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    pentabox = DistevalLibrary('pentabox_fin/disteval/pentabox_fin.json')

    # integrate
    str_result = pentabox(parameters={"s12": 5.0, 's23': -4.0, 's34':2.0, 's45':-6.0, 's51':3.0}, verbose=True, epsrel = 1e-4, epsabs = 1e-100)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')