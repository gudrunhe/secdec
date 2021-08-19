#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    name = "formfactor1L_massive_ebr"
    intlib = IntegralLibrary(f"{name}/{name}_pylink.so")
    intlib.use_Qmc(transform="korobov3", fitfunction="polysingular", verbosity=1)

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = intlib(real_parameters=[100,0.01])

    # convert complex numbers from c++ to sympy notation
    str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')

    # convert result to sympy expressions
    integral_result = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    integral_result_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))

    # examples how to access individual orders
    print('Numerical Result')
    for power in [-2, -1, 0]:
        valreal, valimg = integral_result.coeff('eps',power).coeff('value').as_real_imag()
        errreal, errimg = integral_result.coeff('eps',power).coeff('error').as_real_imag()
        print("eps^{:<2} {: .15f}{:+.15f}*I +/- {:.15e}{:+.15e}*I".format(power,float(valreal),float(valimg),float(errreal),float(errimg)))

    # analytic result
    # from Eq(6.30) of arXiv:1111.2589
    # N[-1/Qsq*(1/2*Log[t]^2 - Log[t]*Log[1 - t] - PolyLog[2, t] + Pi^2/3) /. {Qsq -> 100, t -> 1/100}, 20]
    print('\nAnalytic Result for Qsq = 100, t = msq/Qsq = 0.01')
    print("eps^0  -0.13837355735881397785+0*I")
