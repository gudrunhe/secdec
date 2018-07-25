from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
box2L = IntegralLibrary('box2L_hfkk_a/box2L_hfkk_a_pylink.so')

# choose integrator
box2L.use_Vegas(flags=2,epsrel=10**-5,epsabs=10**-5,maxeval=10**7) # ``flags=2``: verbose --> see Cuba manual

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box2L(real_parameters=[-1.0,-0.8,0.1])

print(str_integral_with_prefactor)
