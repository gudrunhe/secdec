from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
N3 = IntegralLibrary('N3/N3_pylink.so')

# choose integrator
N3.use_Divonne(epsrel=1e-5,epsabs=1e-5,maxeval=10**7,border=1e-8,flags=3)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = N3(together=False)

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

# result from arXiv:1610.07059 (equation 3.2)
print('Reference Result')
print('eps^-2: 1.23370055013617    - 6.20475892887384  * 10^-13 * I')
print('eps^-1: 2.8902545096591976  + 3.875784585038738          * I')
print('eps^0: -0.7785996083247692  - 4.123512600516016          * I')
