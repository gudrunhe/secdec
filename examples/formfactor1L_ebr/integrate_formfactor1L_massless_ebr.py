#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
name = "formfactor1L_massless_ebr"
intlib = IntegralLibrary(f"{name}/{name}_pylink.so")
intlib.use_Qmc(transform="korobov3", fitfunction="polysingular")

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = intlib(real_parameters=[100,1,1,1])

# convert complex numbers from c++ to sympy notation
str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')

# convert result to sympy expressions
integral_result = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
integral_result_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))

# examples how to access individual orders
print('Numerical Result')
for power in [-2, -1, 0]:
    valreal, valimg = integral_result.coeff('eps',power).coeff('value').as_real_imag()
    errreal, errimg = integral_result.coeff('eps',power).coeff('error').as_real_imag()
    print("eps^{:<2} {: .5f}{:+.5f}*I +/- {: .5e}{:+.5e}*I".format(power,float(valreal),float(valimg),float(errreal),float(errimg)))
    
# Expected Result for qsq=100, lsq=1, psq=1 (Compared with Long Chen)
# + ((0,0) +/- (2.93878e-19,0))*eps^-2 + ((-1.6865e-11,8.49456e-11) +/- (1.25602e-07,7.58591e-08))*eps^-1 + ((-0.244975,3.27442e-10) +/- (1.20291e-06,1.07642e-06)) + O(eps)
