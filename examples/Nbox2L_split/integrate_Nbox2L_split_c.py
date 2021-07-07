#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary

if __name__ == "__main__":

    # load c++ library
    box2L = IntegralLibrary('Nbox2L_split_c/Nbox2L_split_c_pylink.so')

    # choose integrator
    box2L.use_Qmc(verbosity=0,minn=15173222401,maxeval=1,transform='korobov3',generatingvectors='cbcpt_cfftw1_6')

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box2L(real_parameters=[-1.0,-0.8,0.1], verbose=True)

    # print results
    print('Numerical Result')
    print('integral without prefactor', str_integral_without_prefactor)
    print('prefactor', str_prefactor)
    print('integral with prefactor', str_integral_with_prefactor)
