from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
box2L_jw = IntegralLibrary('box2L_jw/box2L_jw_pylink.so')

# choose integrator
box2L_jw.use_Vegas(flags=2,epsrel=1e-3,maxeval=10**5)

# integrate non-Euclidean point;
s, t, msq = [2.585, -0.7866, 1.0]

# Note: point requires deformation_parameters_maximum to be reduced (else sign check will fail) 
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box2L_jw([s, t, msq],number_of_presamples=4000,deformation_parameters_maximum=5e-2)

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
print('eps^-4:', integral_with_prefactor.coeff('eps',-4).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-4).coeff('error'), ')')
print('eps^-3:', integral_with_prefactor.coeff('eps',-3).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-3).coeff('error'), ')')
print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

