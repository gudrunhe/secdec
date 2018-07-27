from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
ff = IntegralLibrary('formfactor4L/formfactor4L_pylink.so')

# choose integrator
#ff.use_Vegas(flags=2, maxeval=10**6, mineval=10**6)
ff.use_Qmc(
minn=35*10**5,
minm = 64,

cudablocks=128, cudathreadsperblock=64,
maxnperpackage=64, maxmperpackage=32,

#epsrel=1e-100, epsabs=1e-10, maxeval=10**12,

verbosity=3,

transform='baker',
)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = ff(together=True)

print('Numerical Result')
print(str_integral_without_prefactor)
print(str_prefactor)
print(str_integral_with_prefactor)
