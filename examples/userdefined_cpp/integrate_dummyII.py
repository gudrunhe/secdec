from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
integral = IntegralLibrary('dummyII/dummyII_pylink.so')

# choose integrator
integral.use_Vegas(epsrel=5e-4)

# integrate for a certain value of delta
alpha = 0.77
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral([alpha])

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
print('Numerical Result:')
print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
print('eps^1:', integral_with_prefactor.coeff('eps',1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',1).coeff('error'), ')')
print('\n')

# SecDec-3 result
print('SecDec-3 Result:')
print('eps^-2: 0.0952061338 +- 2.839050321e-06')
print('eps^-1: -2.560694287 +- 0.000569146')
print('eps^0:  21.11963361  +- 0.003467324')
print('eps^1: -78.70007348  +- 0.018435777')
