#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    box = IntegralLibrary('ladder2g/ladder2g_pylink.so')

    # choose integrator
    box.use_Qmc(verbosity=2,minn=10**7,maxeval=1,transform='korobov3')

    # integrator settings used to run the timings
    #box.use_Qmc(verbosity=2,minn=10**6,maxeval=1,transform='korobov3',fitfunction='polysingular')
    #box.use_Qmc(verbosity=2,minn=10**5,maxeval=1,transform='korobov3',fitfunction='polysingular',devices=[-1])
    #box.use_Vegas(flags=2,maxeval=10**6,epsrel=1e-100,epsabs=1e-100)

    # integrate non-Euclidean point;
    # s12, s23, s34, s45, s51, msq, mhsq = [2408664.483714835,-554153.4928655404,1496874.6133096598,735785.915792151,-517267.51076904044,256998.,41487.7]
    s12, s23, s34, s45, s51, msq, mhsq = [512775.5889493924,-75616.7189298124,211810.77151387237,187097.41183131794,-198221.48657428243,256998.,41487.7]
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = \
            box([s12/s12,s23/s12,s34/s12,s45/s12,s51/s12,msq/s12,mhsq/s12],number_of_presamples=10**7)

    # convert complex numbers from c++ to sympy notation
    str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
    str_prefactor = str_prefactor.replace(',','+I*')
    str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')

    # convert result to sympy expressions
    integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    prefactor = sp.sympify(str_prefactor)
    integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
    integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

    # numerical result
    print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
