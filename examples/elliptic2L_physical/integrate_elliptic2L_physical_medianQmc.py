#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    elliptic2L_physical = IntegralLibrary('elliptic2L_physical/elliptic2L_physical_pylink.so')

    # choose integrator
    elliptic2L_physical.use_Qmc(transform='korobov1',verbosity=2, lattice_candidates=11)

    # integrate non-Euclidean point;
    s, t, pp4, msq = [90.,-2.5,1.6,1.]
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = elliptic2L_physical([s, t, pp4, msq], verbose=True)

    # convert complex numbers from c++ to sympy notation
    str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
    str_prefactor = str_prefactor.replace(',','+I*')
    str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')

    # convert result to sympy expressions
    integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    prefactor = sp.sympify(str_prefactor)
    integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
    integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

    # numerical result
    print('Numerical Result for non-Euclidean point:')
    print('does NOT contain prefactor (s/msq)^(3/2)')
    print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

    print('Result from SecDec-3 for non-Euclidean point:')
    print('eps^0: -0.0442887356+I*0.01606818343 +- (2.45625148e-05+I*2.662108194e-05)')

