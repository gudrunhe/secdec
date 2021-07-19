#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
name = "bubble1L_expansion_by_regions_large_mass"
bubble1L = IntegralLibrary(f"{name}/{name}_pylink.so")
bubble1L.use_Qmc(transform="korobov3", fitfunction="polysingular", verbosity=1)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = bubble1L(real_parameters=[0.002,4])

# convert complex numbers from c++ to sympy notation
str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')

# convert result to sympy expressions
integral_result = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
integral_result_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))

# examples how to access individual orders
print('Numerical Result')
for power in [-1, 0]:
    valreal, valimg = integral_result.coeff('eps',power).coeff('value').as_real_imag()
    errreal, errimg = integral_result.coeff('eps',power).coeff('error').as_real_imag()
    print("eps^{:<2} {: .15f}{:+.15f}*I +/- {:.15f}{:+.15f}*I".format(power,float(valreal),float(valimg),float(errreal),float(errimg)))