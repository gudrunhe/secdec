from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
integral = IntegralLibrary('A92/A92_pylink.so')

# choose integrator
integral.use_Divonne(flags=2,epsrel=1e-3,epsabs=1e-10,maxeval=10000000,seed=27, border=1e-8)
#integral.use_Vegas(flags=2,epsrel=1e-3,epsabs=0,nstart=500,nincrease=10000,maxeval=100000000,seed=27,real_complex_together=True) 
#integral.use_Vegas(flags=2,epsrel=1e-5,epsabs=1e-30000,nstart=100000,nincrease=30000,maxeval=10000000) 
#integral.use_Vegas(flags=2,maxeval= 1000000, epsrel=1.e-3, epsabs=0,seed=27)  # fast evaluation
# ``flags=2``: verbose --> see Cuba manual
#integral.use_Cuhre(flags=2,epsrel=1e-2,epsabs=1e-10,maxeval=100000000,real_complex_together=True) 
#integral.use_CQuad(verbose=True,epsrel=1e-3,epsabs=0) 

# integrate
#str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral(number_of_presamples=20000,deformation_parameters_maximum=1e-1,deformation_parameters_minimum=1e-10,together=True)
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral()

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

# Print commands for the full result
print('finite part:', sp.im(integral_with_prefactor.coeff('eps',0).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',0).coeff('error')), ')')
print('1/eps part:', sp.im(integral_with_prefactor.coeff('eps',-1).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',-1).coeff('error')), ')')
print('1/eps^2 part:', sp.im(integral_with_prefactor.coeff('eps',-2).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',-2).coeff('error')), ')')
print('1/eps^3 part:', sp.im(integral_with_prefactor.coeff('eps',-3).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',-3).coeff('error')), ')')
print('1/eps^4 part:', sp.im(integral_with_prefactor.coeff('eps',-4).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',-4).coeff('error')), ')')
print('1/eps^5 part:', sp.im(integral_with_prefactor.coeff('eps',-5).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',-5).coeff('error')), ')')
print('1/eps^6 part:', sp.im(integral_with_prefactor.coeff('eps',-6).coeff('value')), '+/- (', sp.im(integral_with_prefactor_err.coeff('eps',-6).coeff('error')), ')')

