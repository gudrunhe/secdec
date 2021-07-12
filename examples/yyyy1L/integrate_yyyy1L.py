#!/usr/bin/env python3

import math
from pySecDec.integral_interface import IntegralLibrary

# Analytic result (for comparison)
def yyyy1L_analytic(t, u):
    s = -t -u
    ltu = math.log(t/u)
    return -8*( 1 + (t-u)/s*ltu  + (t*t + u*u)/(2*s*s)*\
    ( math.pi*math.pi + ltu**2 ) )

if __name__ == "__main__":

    # load c++ library
    amp = IntegralLibrary('yyyy1L/yyyy1L_pylink.so')

    # choose Qmc integrator
    amp.use_Qmc()

    # integrate
    _, _, result = amp([-1.3,-0.8]) # t, u

    # print result
    print('Numerical Result:', result)
    print('Analytic Result:', yyyy1L_analytic(-1.3,-0.8))
