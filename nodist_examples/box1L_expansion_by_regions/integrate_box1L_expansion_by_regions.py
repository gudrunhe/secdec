from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
name = "box1L_expansion_by_regions"
intlib = IntegralLibrary("{0}/{0}_pylink.so".format(name))
intlib.use_Qmc(transform="korobov3", fitfunction="polysingular")

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = intlib(real_parameters=[4.0, -2.82842712475, -2.82842712475, 0.1])

# convert complex numbers from c++ to sympy notation
str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')

# convert result to sympy expressions
integral_result = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
integral_result_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))

# examples how to access individual orders
print('Numerical Result')
for power in [0]:
    valreal, valimg = integral_result.coeff('eps',power).coeff('value').as_real_imag()
    errreal, errimg = integral_result.coeff('eps',power).coeff('error').as_real_imag()
    print("eps^{:<2} {: .5f}{:+.5f}*I +/- {: .5f}{:+.5f}*I".format(power,float(valreal),float(valimg),float(errreal),float(errimg)))
