from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
elliptic2L_euclidean = IntegralLibrary('elliptic2L_euclidean/elliptic2L_euclidean_pylink.so')

# choose integrator
elliptic2L_euclidean.use_Vegas(epsrel=1e-5,maxeval=10**7)

# integrate
s, t, pp4, msq = [-4./3.,-16./5.,-100./39.,1.] # Euclidean point
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = elliptic2L_euclidean([s, t, pp4, msq])

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
print('Numerical Result for Euclidean point')
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

# elliptic Integral f^A_66 of Eq(4.21) arXiv:1609.06685
# analytic result kindly provided by the authors
print('Analytic Result for Euclidean point')
print('eps^0: 0.247074199140732131068066')
