#!/usr/bin/env python3

from pySecDec.integral_interface import DistevalLibrary
from math import log

if __name__ == "__main__":

    # load c++ library
    easy = DistevalLibrary('easy/disteval/easy.json')

    # integrate
    result = easy(epsrel=1e-5)

    # print result
    print('Numerical Result:', result)
    print('Analytic Result:' + ' + (%f)*eps^-1 + (%f) + O(eps)' % (1.0,1.0-log(2.0)))
