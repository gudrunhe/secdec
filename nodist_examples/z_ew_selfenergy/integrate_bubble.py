from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import os

if __name__ == "__main__":

    # load c++ library
    integral = IntegralLibrary('bubble/bubble_pylink.so')

    # choose integrator
    integral.use_Qmc(verbosity=2,minn=10**5,maxeval=0, minm=10 ,transform='korobov3')

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral(real_parameters=[-0.01],
    number_of_presamples=10000, deformation_parameters_maximum=0.01, deformation_parameters_minimum=0.00000001)

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

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
    print('eps^1 :', integral_with_prefactor.coeff('eps',1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',1).coeff('error'), ')')
