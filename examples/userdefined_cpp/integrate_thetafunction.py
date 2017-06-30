from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
integral = IntegralLibrary('thetafunction/thetafunction_pylink.so')

# choose integrator
integral.use_Vegas(epsrel=1e-4, maxeval=10**7)

# integrate for a certain value of delta
delta = 0.1
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral([delta])

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
print('Numerical Result')
print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
print('eps^1:', integral_with_prefactor.coeff('eps',1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',1).coeff('error'), ')')

# mathematica result
print('Mathematica Result')
print('eps^-2: 0')
print('eps^-1: -0.42618702305960626')
print('eps^0: -0.5909427273863541')
print('eps^1: -0.7974333289506665')
