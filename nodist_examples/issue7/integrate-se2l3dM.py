#! /usr/bin/env python3
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    se2l3dM = IntegralLibrary('se2l3dM/se2l3dM_pylink.so')
    #se2l3dM.use_Cuhre(flags=0, epsrel=1e-4)

    m=0
    for s in [-0.1, -0.5, -1.0, -1.5, -2.0]:
        print('J({s}) = {res}'.format(s=s, res=se2l3dM([s,m])[2]))

