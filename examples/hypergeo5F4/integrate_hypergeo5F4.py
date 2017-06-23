from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
hypergeo5F4 = IntegralLibrary('hypergeo5F4/hypergeo5F4_pylink.so')

# choose integrator
hypergeo5F4.use_Vegas(flags=2, epsrel=1e-4) # ``flags=2``: verbose --> see Cuba manual

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = hypergeo5F4([0.5])

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
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
print('eps^1:', integral_with_prefactor.coeff('eps',1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',1).coeff('error'), ')')
print('eps^2:', integral_with_prefactor.coeff('eps',2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',2).coeff('error'), ')')
print('eps^3:', integral_with_prefactor.coeff('eps',3).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',3).coeff('error'), ')')
print('eps^4:', integral_with_prefactor.coeff('eps',4).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',4).coeff('error'), ')')

# analytic result from HypExp (D.Maitre, T.Huber, hep-ph/0507094)
print('Analytic Result')
print('eps^0: 1')
print('eps^1: 0.1895324')
print('eps^2: -2.2990427')
print('eps^3: 55.469019')
print('eps^4: -1014.3924')
