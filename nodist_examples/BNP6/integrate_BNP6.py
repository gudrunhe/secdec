#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
from pySecDec.integral_interface import series_to_ginac, series_to_sympy, series_to_mathematica, series_to_maple

import sympy as sp

if __name__ == "__main__":

    # load c++ library
    BNP6 = IntegralLibrary('BNP6/BNP6_pylink.so')

    # choose integrator
    BNP6.use_Qmc()

    # integrate
    _, _, str_integral_with_prefactor = BNP6(real_parameters=[9.,-2.5],verbose=True)
    integral_with_prefactor, integral_with_prefactor_err = series_to_sympy(str_integral_with_prefactor)
    integral_with_prefactor = sp.sympify(integral_with_prefactor)
    integral_with_prefactor_err = sp.sympify(integral_with_prefactor_err)

    # print result
    print('Numerical Result')
    print('eps^-2 with prefactor:', integral_with_prefactor.coeff('eps',-2), '+/- (', integral_with_prefactor_err.coeff('eps',-2),')')
    print('eps^-1 with prefactor:', integral_with_prefactor.coeff('eps',-1), '+/- (', integral_with_prefactor_err.coeff('eps',-1), ')')
    print('eps^0 with prefactor:', integral_with_prefactor.coeff('eps',0), '+/- (', integral_with_prefactor_err.coeff('eps',0), ')')

    print('Analytic Result at (s,t,u)=(9,-2.5,-6.5)')
    print('eps^-2 with prefactor:', '0.10907856854318447 - 0.5799863360473464*I')
    print('eps^-1 with prefactor:', '-0.8876663743916553 + 4.360251717854891*I')
    print('eps^0 with prefactor: ', '0.7966721383373115 - 18.22048104236002*I')
