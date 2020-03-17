from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
BNP6 = IntegralLibrary('BNP6/BNP6_pylink.so')

# choose integrator
#BNP6.use_Qmc(verbosity=3,devices=[-1,0,1,2,3],minn=10**8,transform='korobov3') 
BNP6.use_Vegas(flags=2,epsrel=1e-2,epsabs=1e-10,nstart=10000,nincrease=1000,maxeval=10000000) # 

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = BNP6(real_parameters=[9.,-2.5,-6.5])

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

# print result
print('Numerical Result')
print('prefactor', prefactor)
#print('eps^0 without prefactor:', integral_without_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_without_prefactor_err.coeff('eps',0).coeff('error'), ')')
print('eps^-2 with prefactor:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1 with prefactor:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0 with prefactor:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

print('Analytic Result at (s,t,u)=(9,-2.5,-6.5)')
print('eps^-2 with prefactor:', '0.10907856854318447 - 0.5799863360473464*I')
print('eps^-1 with prefactor:', '-0.8876663743916553 + 4.360251717854891*I')
print('eps^0 with prefactor: ', '0.7966721383373115 - 18.22048104236002*I')
