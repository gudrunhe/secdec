from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
box2L = IntegralLibrary('box2L/box2L_pylink.so')

# choose integrator
box2L.use_Vegas(flags=2,epsrel=1e-3,maxeval=10**7,nstart=1000,nincrease=500) # Gives wrong result (!)
#box2L.use_Vegas(flags=2,epsrel=1e-3,maxeval=10**7,nstart=10000,nincrease=5000) # Gives correct result
#box2L.use_Qmc(verbosity=3,minn=100000,maxeval=1,transform='korobov3') # Gives very large error

# Note: point evaluates fine with default settings
#str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box2L([100.0, -29.17011048, 2.473471074]) # s=100000000, t=250628., mt=173

# Note: point requires nstart > 1000
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box2L([100.0, -0.2506281447, 0.02992900000]) # s=100000000, t=250628., mt=173

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
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

