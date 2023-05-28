#!/usr/bin/env python3

from pySecDec.integral_interface import DistevalLibrary, IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    disteval_lib = DistevalLibrary('hh_b110/disteval/hh_b110.json')
    integral_lib = IntegralLibrary('hh_b110/hh_b110_pylink.so')

    s12, s13, mt2 = [3.1,1.8,1]
    disteval_result = disteval_lib(parameters={'s12': s12, 's13': s13, 'mt2' : mt2})
    _, _, integral_result = integral_lib(real_parameters=[ s12, s13, mt2])
   
    print(disteval_result)
    print(integral_result)

