from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary

# load c++ library
bubble_6L = IntegralLibrary('bubble6L/bubble6L_pylink.so')

# choose integrator
bubble_6L.use_Qmc(verbosity=2,minn=10**7,maxeval=1,transform='korobov3',fitfunction='polysingular')

# integrator settings used to run the timings
#bubble_6L.use_Qmc(verbosity=2,minn=10**6,maxeval=1,transform='baker',fitfunction='polysingular')
#bubble_6L.use_Qmc(verbosity=2,minn=4*10**4,maxeval=1,transform='baker',fitfunction='polysingular',devices=[-1])
#bubble_6L.use_Vegas(flags=2,maxeval=10**5,epsrel=1e-100,epsabs=1e-100)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = bubble_6L()

# examples how to access individual orders
print('Numerical Result')
print('integral without prefactor', str_integral_without_prefactor)
print('prefactor', str_prefactor)
print('integral with prefactor', str_integral_with_prefactor)

# analytic result from arXiv:1705.06483
# 9.264208985946416/eps^2 + 91.73175282208716/eps + O(eps^0)
