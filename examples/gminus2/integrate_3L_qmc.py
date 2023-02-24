#!/usr/bin/env python3
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    vertex3L = IntegralLibrary('g2_3L/g2_3L_pylink.so')

    # choose integrator
    vertex3L.use_Qmc()

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = vertex3L(real_parameters=[1.0], verbose=True)

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
 #   print('Numerical Result')
#    print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

    print('Analytic Result')
    print('eps^-2: -0.1428571429')
    print('eps^-1: 0.6384337090')
    print('eps^0 : -0.426354612+I*1.866502363')
