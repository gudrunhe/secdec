from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary

# load c++ library
banana_3mass = IntegralLibrary('banana_3mass/banana_3mass_pylink.so')

# choose integrator
banana_3mass.use_Qmc(verbosity=3,minn=10**6,maxeval=10**6,transform='korobov2')

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = banana_3mass(real_parameters=[30.0,1.0,1.3,0.7], deformation_parameters_maximum=0.1)

# print results
print('Numerical Result')
print('integral without prefactor', str_integral_without_prefactor)
print('prefactor', str_prefactor)
print('integral with prefactor', str_integral_with_prefactor)

