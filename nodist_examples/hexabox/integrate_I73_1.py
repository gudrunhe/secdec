#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary, series_to_sympy
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    I73_1 = IntegralLibrary('I73_1/I73_1_pylink.so')

    # choose integrator
    #I73_1.use_Qmc(verbosity=3,devices=[-1,0,1,2,3],minn=10**8,transform='korobov3')
    I73_1.use_Vegas(flags=2,epsrel=1e-3,maxeval=10**7)


    # integrate non-Euclidean point;
    v1,v2,v3,v4,v5 = [-3.,-3.,-1.,-1.,-1.]
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = I73_1([v1,v2,v3,v4,v5])

    # convert result to sympy expressions
    integral_with_prefactor, integral_with_prefactor_err = series_to_sympy(str_integral_with_prefactor)
    integral_with_prefactor = sp.sympify(integral_with_prefactor)
    integral_with_prefactor_err = sp.sympify(integral_with_prefactor_err)

    # numerical result
    print('eps^0:', integral_with_prefactor.coeff('eps',0), '+/- (', integral_with_prefactor_err.coeff('eps',0), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1), '+/- (', integral_with_prefactor_err.coeff('eps',-1), ')')
    print('eps^-2:', integral_with_prefactor.coeff('eps',-2), '+/- (', integral_with_prefactor_err.coeff('eps',-2), ')')
    print('eps^-3:', integral_with_prefactor.coeff('eps',-3), '+/- (', integral_with_prefactor_err.coeff('eps',-3), ')')
    print('eps^-4:', integral_with_prefactor.coeff('eps',-4), '+/- (', integral_with_prefactor_err.coeff('eps',-4), ')')
