#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    box1L = DistevalLibrary('box1L/disteval/box1L.json')

    # integrate
    result = box1L(parameters={"s": 4.0, "t": -0.75, "s1": 1.25, "msq": 1.0}, verbose=False, format="json")
    values = result["sums"]["box1L"]

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', values[(-2,)][0], '+/- (', values[(-2,)][1], ')')
    print('eps^-1:', values[(-1,)][0], '+/- (', values[(-1,)][1], ')')
    print('eps^0 :', values[( 0,)][0], '+/- (', values[( 0,)][1], ')')

    print('Analytic Result')
    print('eps^-2: -0.1428571429')
    print('eps^-1: 0.6384337090')
    print('eps^0 : -0.426354612+I*1.866502363')
