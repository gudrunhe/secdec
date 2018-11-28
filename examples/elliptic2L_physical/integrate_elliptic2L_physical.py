from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
elliptic2L_physical = IntegralLibrary('elliptic2L_physical/elliptic2L_physical_pylink.so')

# choose integrator
elliptic2L_physical.use_Qmc(verbosity=2,minn=2147483647,maxeval=1,transform='korobov1',fitfunction='polysingular')

# integrator settings used to run the timings
#elliptic2L_physical.use_Qmc(verbosity=2,minn=3*10**6,maxeval=3*10**6,transform='korobov1',fitfunction='polysingular')
#elliptic2L_physical.use_Qmc(verbosity=2,minn=3*10**6,maxeval=3*10**6,transform='korobov1',fitfunction='polysingular',devices=[-1])
#elliptic2L_physical.use_Vegas(flags=2,maxeval=3*10**7,epsrel=1e-100,epsabs=1e-100)
 
# integrate non-Euclidean point;
s, t, pp4, msq = [90.,-2.5,1.6,1.]
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = elliptic2L_physical([s, t, pp4, msq])

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
print('Numerical Result for non-Euclidean point:')
print('does NOT contain prefactor (s/msq)^(3/2)')
print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

print('Result from SecDec-3 for non-Euclidean point:')
print('eps^0: -0.0442887356+I*0.01606818343 +- (2.45625148e-05+I*2.662108194e-05)')

