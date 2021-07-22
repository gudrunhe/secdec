#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
from math import log

if __name__ == "__main__":

    #fix expansion parameter
    t = 0.01

    # load c++ library
    intlib = IntegralLibrary('make_regions/make_regions_pylink.so')

    # integrate
    _, _, result = intlib(real_parameters=[t])

    # print result
    print('Numerical Result:' + result)
    print('Analytic Result:' + ' (%f) + O(delta)' % (-log(t)))
