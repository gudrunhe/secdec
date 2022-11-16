#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    box1L = DistevalLibrary('box1L/disteval/box1L.json')

    # integrate
    str_result = box1L(parameters={"s": 4.0, "t": -0.75, "s1": 1.25, "msq": 1.0}, verbose=False)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')

    print('Analytic Result')
    print('eps^-2: -0.1428571429')
    print('eps^-1: 0.6384337090')
    print('eps^0 : -0.426354612+I*1.866502363')
