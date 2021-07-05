#!/usr/bin/env python3

import pySecDec.integral_interface as ii
import sympy as sp

# load c++ library
d0_direct = ii.IntegralLibrary(f"D0_direct/D0_direct_pylink.so")
d0_reduced = ii.IntegralLibrary(f"D0_reduced/D0_reduced_pylink.so")

d0_direct.use_Qmc(transform="korobov3", fitfunction="polysingular")
d0_reduced.use_Qmc(transform="korobov3", fitfunction="polysingular")

def pp(integral_value, powers):
    value, stdev = ii.series_to_sympy(str_integral)
    value_sp = sp.sympify(value)
    stdev_sp = sp.sympify(stdev)
    for power in powers:
        v_re, v_im = value_sp.coeff('eps',power).as_real_imag()
        s_re, s_im = stdev_sp.coeff('eps',power).as_real_imag()
        print(f"eps^{power:+d}: " \
            f"{float(v_re):22.15e}{float(v_im):+22.15e}*I +/- " \
            f"{float(s_re):8.2e}{float(s_im):+8.2e}*I")

# integrate
print("Direct:")
str_integral_without_prefactor, str_prefactor, str_integral = d0_direct(real_parameters=[1.0, -0.4])
pp(str_integral, [-2,-1,0])

print("Reduced:")
str_integral_without_prefactor, str_prefactor, str_integral = d0_reduced(real_parameters=[1.0, -0.4])
pp(str_integral, [-2,-1,0])
