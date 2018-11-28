from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary

# load c++ library
ff = IntegralLibrary('formfactor4L/formfactor4L_pylink.so')

# choose integrator
ff.use_Qmc(
minn=35*10**5,

# settings used for timings
#minn=5*10**3, devices = [-1],
#minn=2*10**5,

maxeval=1,
minm = 64,

cudablocks=128, cudathreadsperblock=64,
maxnperpackage=8, maxmperpackage=8,

verbosity=2,

transform='baker',

fitfunction='polysingular',
)
#ff.use_Vegas(flags=2,maxeval=10**5,epsrel=1e-100,epsabs=1e-100) # used for timings

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = ff()

# print results
print('Numerical Result')
print('integral without prefactor', str_integral_without_prefactor)
print('prefactor', str_prefactor)
print('integral with prefactor', str_integral_with_prefactor)
