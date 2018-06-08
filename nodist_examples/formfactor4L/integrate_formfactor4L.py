from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
ff = IntegralLibrary('formfactor4L/formfactor4L_pylink.so')

# choose integrator
#ff.use_Vegas(flags=2, maxeval=10**6, mineval=10**6)
ff.use_Qmc(
minn=15*10**7,
minm = 64,

cudablocks=512, cudathreadsperblock=256,
maxnperpackage=64, maxmperpackage=2**10,

epsrel=1e-100, epsabs=1e-10, maxeval=10**12,

verbosity=2,
)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = ff(together=False)

print('Numerical Result')
print(str_integral_without_prefactor)
print(str_prefactor)
print(str_integral_with_prefactor)
