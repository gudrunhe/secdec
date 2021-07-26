#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    bubble2L_full = IntegralLibrary('bubble2L_full/bubble2L_full_pylink.so')

    # choose integrator
    bubble2L_full.use_Vegas(flags=2,epsrel=1e-3,epsabs=1e-10,nstart=5000, nincrease=10000, maxeval=10000000) # ``flags=2``: verbose --> see Cuba manual

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = bubble2L_full(number_of_presamples=1000000,deformation_parameters_maximum=0.5,real_parameters=[0.003, 1.0])

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

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

    # analytic result as expansion: only converges for s<<msq !!!

    #I3[s_,msq_,nmax_]:= Module[{splusidelta,f0,f1},
    #splusidelta=s+I*10^(-15);
    #f0 = 3/2/msq*Sum[Factorial[n]*Gamma[3/2]/(Gamma[n+3/2]*(n+1)^2)*(-s/4/msq)^n,{n,0,nmax}];
    #f1 = -1/2/msq*Sum[Factorial[n]*Gamma[3/2]/(Gamma[n+3/2]*(n+1))*(-s/4/msq)^n,{n,0,nmax}];
    #res =f0+f1*Log[-splusidelta/msq];
    #Return[res]
    #]

    # for s=0.003, msq=1:
    # mathematica result,nmax=1:   4.403657852283224 + 1.570403627713031*I
    # mathematica result,nmax=2:   4.403658192740373 + 1.5704037847926637*I
    # mathematica result,nmax=100: 4.403658192582334 + 1.5704037847169694*I
    # pySecDec result full:        4.40545611430997308 + 1.57124257362624387*I +/- ( 0.00432238468534031951 + 0.00153966152208577533*I )
    # pySecDec expansion order 1:  4.40544,1.57117) +/- (0.000367815,0.000152761) + ((-5.45655e-11,-1.07502e-11) +/- (6.75993e-05,3.05832e-06))*eps^-1
