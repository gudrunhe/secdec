from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary

# load c++ library
box2L = IntegralLibrary('Nbox2L_split_a/Nbox2L_split_a_pylink.so')

# choose integrator
box2L.use_Qmc(verbosity=2,minn=10**7,maxeval=1,transform='korobov4',fitfunction='polysingular')

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box2L(real_parameters=[-1.0,-0.8,0.1])

# print results
print('Numerical Result')
print('integral without prefactor', str_integral_without_prefactor)
print('prefactor', str_prefactor)
print('integral with prefactor', str_integral_with_prefactor)
