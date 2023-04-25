#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary

lib = IntegralLibrary('b25_8p8d-cpu/b25_8p8d_pylink.so')

lib.use_Qmc(verbosity=0, minn=10000)

p = [12/23, 1, 262/35, -145/53, 327/164, -101/36, 249/124]
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = lib(p, epsrel=1e-3, epsabs=1e-10, verbose=True, number_of_presamples=10000)

print(str_integral_with_prefactor)
