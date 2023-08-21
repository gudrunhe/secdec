#!/usr/bin/env python3
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
box1L = IntegralLibrary('box1L/box1L_pylink.so')

# choose integrator
box1L.use_Qmc()

# integrate
result_without_prefactor, result_prefactor, result_with_prefactor = \
    box1L(real_parameters=[4.0, -0.75, 1.25, 1.0],
          epsrel=1e-3, epsabs=1e-10, format="json")
values = result_with_prefactor["sums"]["box1L"]

# examples how to access individual orders
print('Numerical Result')
print('eps^-2:', values[(-2,)][0], '+/- (', values[(-2,)][1], ')')
print('eps^-1:', values[(-1,)][0], '+/- (', values[(-1,)][1], ')')
print('eps^0 :', values[( 0,)][0], '+/- (', values[( 0,)][1], ')')

print('Analytic Result')
print('eps^-2: -0.1428571429')
print('eps^-1: 0.6384337090')
print('eps^0 : -0.426354612+I*1.866502363')
