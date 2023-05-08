#!/usr/bin/env python3

from pySecDec.integral_interface import DistevalLibrary

lib = DistevalLibrary('hexatriangle/disteval/hexatriangle.json')

p = {"mh2": 12/23, 'mt2': 1, 'x12': 262/35, 'x23': -145/53, 'x35': 327/164, 'x41': -101/36, 'x54': 249/124}
result = lib(p, points=10000, number_of_presamples=10000, epsrel=1e-3, epsabs=1e-10, verbose=True)

print("Known value: 1.454919812(7)*10^-7 - 1.069797219(8)*10^-07 I + O(eps)")
print("Calculated value:", result)
