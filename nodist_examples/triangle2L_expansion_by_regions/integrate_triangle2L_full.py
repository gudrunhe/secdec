from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
triangle2L_case8_full = IntegralLibrary('triangle2L_case8_full/triangle2L_case8_full_pylink.so')

# choose integrator
triangle2L_case8_full.use_Vegas(flags=2,epsrel=1e-3,epsabs=1e-10,nstart=5000, nincrease=10000, maxeval=10000000) # ``flags=2``: verbose --> see Cuba manual

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = triangle2L_case8_full(number_of_presamples=1000000,deformation_parameters_maximum=0.5,real_parameters=[2.0, 1.0])

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
print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

#print('Analytic Result')
#print('eps^-2: ')
#print('eps^-1: ')
#print('eps^0 : ')

# double pole: (Mathematica syntax)

#P2[s_,msq_]:= Module[{splusidelta},
#splusidelta=s+I*10^(-15);
#res=1/s^2 *(Log[-splusidelta/msq]*Log[1+splusidelta/msq]+PolyLog[2,-splusidelta/msq]);
#Return[res]
#]

