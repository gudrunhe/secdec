#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary

lib = IntegralLibrary('b25_8p8d-cpu/b25_8p8d_pylink.so')

lib.use_Qmc(verbosity=0, minn=10000)

p = [12/23, 1, 262/35, -145/53, 327/164, -101/36, 249/124]
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = lib(p, epsrel=1e-3, epsabs=1e-10, verbose=True, number_of_presamples=10000)

print("Known value: 1.454919812(7)*10^-7 - 1.069797219(8)*10^-07 I + O(eps)")
print("Calculated value:", str_integral_with_prefactor)
