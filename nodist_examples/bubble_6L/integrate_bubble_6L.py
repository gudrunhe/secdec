from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
bubble_6L = IntegralLibrary('bubble_6L/bubble_6L_pylink.so')

# choose integrator
#bubble_6L.use_Qmc(verbosity=3,devices=[-1,0,1,2,3],minn=10**8) 
bubble_6L.use_Vegas(flags=2,epsrel=1e-3,epsabs=1e-8,nstart=10000,nincrease=1000,maxeval=10000000) # 

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = bubble_6L(real_parameters=[-1.])

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
print('prefactor', prefactor)
print('eps^-2 with prefactor:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
print('eps^-1 with prefactor:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
print('eps^0 with prefactor:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

# analytic result: 

# epsm2 = 147/16*Zeta[7];
# epsm1 = -(147/16*Zeta[7]+27/2*Zeta[3]*Zeta[5]+27/10*zeta35-2063/504000*Pi^8);
# eps0 = unknown
# zeta35=0.037707673;
# prefac=Gamma[1-eps]^2*Gamma[1+eps]/Gamma[2-2*eps];

# pysecdec prefactor (from Feynman parametrisation) is Gamma[6*eps]
# prefactor in 1705.06483 eq.(A3) is prefac=Gamma[1-eps]^2*Gamma[1+eps]/Gamma[2-2*eps] per loop

# result = N[Normal[Series[prefac^6*(epsm2/eps^2+epsm1/eps+fin), {eps,0,0}]]];
# result = 9.264208985946416/eps^2 + 91.73175282208716/eps + finite

# pysecdec result:
# eps^-2: 9.26430342193453171 +/- ( 0.000899765649050767845 )
# eps^-1: 91.7338323848475738 +/- ( 0.0125856346079082672 )
# eps^0: 1118.67438974297716 +/- ( 0.132744000721334243 )
