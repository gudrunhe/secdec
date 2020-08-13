from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    I73_29 = IntegralLibrary('I73_29/I73_29_pylink.so')

    # choose integrator
    #I73_29.use_Qmc(verbosity=3,devices=[-1,0,1,2,3],minn=10**8,transform='korobov3')
    I73_29.use_Vegas(flags=2,epsrel=1e-3,maxeval=10**7)


    # integrate non-Euclidean point;
    v1,v2,v3,v4,v5 = [-3.,-3.,-1.,-1.,-1.]
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = I73_29([v1,v2,v3,v4,v5])

    # convert result to sympy expressions
    integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    prefactor = sp.sympify(str_prefactor)
    integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
    integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

    # numerical result
    print('eps^0:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('eps^-3:', integral_with_prefactor.coeff('eps',-3).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-3).coeff('error'), ')')
    print('eps^-4:', integral_with_prefactor.coeff('eps',-4).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-4).coeff('error'), ')')
