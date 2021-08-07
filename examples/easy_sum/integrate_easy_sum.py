#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
from math import log

if __name__ == "__main__":

    # load c++ library
    easy = IntegralLibrary('easy_sum/easy_sum_pylink.so')

    # integrate
    _, _, result = easy([0.1]) # s

    # print result
    print('Numerical Result:' + result)
    print('Analytic Result:')
    print('+ 0.2/eps + 0.22962348064032504712')
    print('0.05/eps^2 + 0.015342640972002734529/eps + 0.0033313156240476989125')
