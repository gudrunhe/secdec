from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
triangle2L = IntegralLibrary('ff3_massive/ff3_massive_pylink.so')

# choose integrator
#triangle2L.use_Vegas(nstart=10**4, maxeval=10**8, flags=2, real_complex_together=True) # ``flags=2``: verbose --> see Cuba manual
#triangle2L.use_Suave(maxeval=10**5, flags=2, real_complex_together=True) # ``flags=2``: verbose --> see Cuba manual
triangle2L.use_Qmc(verbosity=3, cudablocks=1024, cudathreadsperblock=256, minn=10**7, transform='baker',epsrel=1e-4, epsabs=1e-7,maxeval=10**6,maxmperpackage=16)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = triangle2L(real_parameters=[0.72135,1.0],number_of_presamples=1000000,deformation_parameters_maximum=0.1,deformation_parameters_minimum=1e-10, together=True)

print('Numerical Result')
print(str_integral_without_prefactor)
print(str_prefactor)
print(str_integral_with_prefactor)
