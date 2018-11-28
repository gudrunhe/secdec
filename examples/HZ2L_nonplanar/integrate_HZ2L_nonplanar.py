from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
HZ2L_nonplanar = IntegralLibrary('HZ2L_nonplanar/HZ2L_nonplanar_pylink.so')

# choose integrator
HZ2L_nonplanar.use_Qmc(verbosity=2,minn=10**8,maxeval=1,transform='korobov3')

# integrator settings used to run the timings
#HZ2L_nonplanar.use_Qmc(verbosity=2,minn=10**6,maxeval=10**6,transform='korobov3')
#HZ2L_nonplanar.use_Qmc(verbosity=2,minn=10**5,maxeval=10**5,transform='korobov3',devices=[-1])
#HZ2L_nonplanar.use_Vegas(flags=2,maxeval=5*10**5,epsrel=1e-100,epsabs=1e-100)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = HZ2L_nonplanar(real_parameters=[200.,-23.,9.,1.56,0.81])

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
print('eps^-2 with prefactor:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1 with prefactor:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0 with prefactor:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
