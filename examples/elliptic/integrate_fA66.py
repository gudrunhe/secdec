from __future__ import print_function
from __future__ import division
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
fA66 = IntegralLibrary('fA66/fA66_pylink.so')

# choose integrator
fA66.use_Vegas(flags=2,epsrel=1e-4) # ``flags=2``: verbose --> see Cuba manual

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = fA66(real_parameters=[-4./3.,-16./5.,-100./39.,1.])

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
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

