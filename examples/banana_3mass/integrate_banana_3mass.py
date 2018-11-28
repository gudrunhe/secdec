from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary

# load c++ library
banana_3mass = IntegralLibrary('banana_3mass/banana_3mass_pylink.so')

# choose integrator
banana_3mass.use_Qmc(verbosity=2,minn=10**6,maxeval=1,transform='korobov2')

# integrator settings used to run the timings
#banana_3mass.use_Qmc(verbosity=2,minn=10**5,maxeval=10**5,transform='korobov2')
#banana_3mass.use_Qmc(verbosity=2,minn=10**5,maxeval=10**5,transform='korobov2',devices=[-1])
#banana_3mass.use_Vegas(flags=2,maxeval=10**6,epsrel=1e-100,epsabs=1e-100)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = banana_3mass(real_parameters=[30.0,1.0,1.3,0.7], deformation_parameters_maximum=0.1)

# print results
print('Numerical Result')
print('integral without prefactor', str_integral_without_prefactor)
print('prefactor', str_prefactor)
print('integral with prefactor', str_integral_with_prefactor)

