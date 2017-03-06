from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
formfactor = IntegralLibrary('formfactor3L/formfactor3L_pylink.so')

# choose integrator
formfactor.use_Vegas(epsrel=1e-3)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = formfactor()

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

# analytic result
# 1510.06758 (von Manteuffel, Panzer, Schabinger) anc/MIexpansions.m
# INT["B3", 7, 1722, 7, 0, {0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0}]
print('Analytic Result')
print('eps^0: -34.096929755001252682')
print('eps^1: -295.87002604581508543')
print('eps^2: -2052.9323782881109903')
print('eps^3: -10599.637376067276383')
print('eps^4: -49873.467692019361424')
