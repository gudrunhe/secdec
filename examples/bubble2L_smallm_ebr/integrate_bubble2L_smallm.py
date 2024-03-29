#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary, series_to_sympy
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    name = "bubble2L_smallm"
    intlib = IntegralLibrary(f"{name}/{name}_pylink.so")
    intlib.use_Qmc(transform="korobov3", verbosity=1)

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = intlib(real_parameters=[4,0.002])
    
    print(str_integral_with_prefactor)

    # convert the result to sympy expressions
    result, error = map(sp.sympify, series_to_sympy(str_integral_with_prefactor))

    # examples how to access individual orders
    print('Numerical Result')
    for power in [-1, 0]:
        valreal, valimg = result.coeff('eps',power).as_real_imag()
        errreal, errimg = error.coeff('eps',power).as_real_imag()
        print("eps^{:<2} {: .15f}{:+.15f}*I +/- {:.15f}{:+.15f}*I".format(power,float(valreal),float(valimg),float(errreal),float(errimg)))
