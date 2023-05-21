#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    g2_4L = DistevalLibrary('g2_4L/disteval/g2_4L.json')

    # integrate
    str_result = g2_4L(parameters={"msq": 1.0}, verbose=True)

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-4:', value.coeff('eps',-4), '+/- (', error.coeff('eps',-4), ')')
    print('eps^-3:', value.coeff('eps',-3), '+/- (', error.coeff('eps',-3), ')')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')
