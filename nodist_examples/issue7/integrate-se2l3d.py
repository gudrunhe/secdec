#! /usr/bin/env python3
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    se2l3d = IntegralLibrary('se2l3d/se2l3d_pylink.so')
    #se2l3d.use_Cuhre(flags=0, epsrel=1e-4)

    for s in [-0.1, -0.5, -1.0, -1.5, -2.0]:
        print('J({s}) = {res}'.format(s=s, res=se2l3d([s])[2]))

