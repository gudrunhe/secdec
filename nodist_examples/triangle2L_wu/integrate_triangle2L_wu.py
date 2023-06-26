from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
from pySecDec.integral_interface import series_to_ginac, series_to_sympy, series_to_mathematica, series_to_maple

if __name__ == "__main__":

    # load c++ library
    triangle2L_wu = IntegralLibrary('triangle2L_wu/triangle2L_wu_pylink.so')

    # choose integrator
    triangle2L_wu.use_Qmc(verbosity=2,maxeval=10**8,epsrel=0.001,epsabs=10**(-7),transform='korobov3')

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = triangle2L_wu()

 
    print(str_integral_with_prefactor)
    print(series_to_sympy(str_integral_with_prefactor))

    # result from arXiv:1610.07059 (equation 3.2)
    # Note: The result given in the reference above has a sign error in the finite part.
    #       The result given below has been confirmed by the authors of arXiv:1610.07059
    #       in a private communication.
    print('Reference Result')
    print('eps^-2: 1.23370055013617    - 6.20475892887384  * 10^-13 * I')
    print('eps^-1: 2.8902545096591976  + 3.875784585038738          * I')
    print('eps^0:  0.7785996083247692  + 4.123512600516016          * I')
