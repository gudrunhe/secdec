from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
from math import pi
import sympy as sp

# load c++ library
two_regulators = IntegralLibrary('two_regulators/two_regulators_pylink.so')

# choose integrator
two_regulators.use_Vegas(flags=2) # ``flags=2``: verbose --> see Cuba manual

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = two_regulators()

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
print('eps^-1 * alpha^-1:', integral_with_prefactor.coeff('alpha',-1).coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('alpha',-1).coeff('eps',-1).coeff('error'), ')')
print('eps^-2 * alpha^0:', integral_with_prefactor.coeff('alpha',0).coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('alpha',-2).coeff('eps',0).coeff('error'), ')')
print('eps^0 * alpha^0:', integral_with_prefactor.coeff('alpha',0).coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('alpha',0).coeff('eps',0).coeff('error'), ')')

# analytic result
print('Analytic Result')
print('eps^-1 * alpha^-1: -1')
print('eps^-2 * alpha^0: 0.5')
print('eps^0 * alpha^0: ' + str(-pi**2/6) )

# we thank Guido Bell and Rudi Rahn for this example
