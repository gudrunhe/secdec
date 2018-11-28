from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
hyperelliptic = IntegralLibrary('hyperelliptic/hyperelliptic_pylink.so')

# choose integrator
hyperelliptic.use_Qmc(verbosity=2,maxeval=1,minn=10**8,transform='korobov3',fitfunction='polysingular')

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = hyperelliptic(real_parameters=[10.0, -0.75, 1.0, 1.3, 0.7])

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
print('prefactor', prefactor)
print('eps^0 without prefactor:', integral_without_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_without_prefactor_err.coeff('eps',0).coeff('error'), ')')
print('eps^0 with prefactor:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

