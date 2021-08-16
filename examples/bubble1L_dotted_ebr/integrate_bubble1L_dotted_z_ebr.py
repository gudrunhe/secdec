#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary, series_to_sympy
import sympy as sp
import numpy as np

# the exact analytical result of the diagram
def f(psq,msq):
    msq = complex(msq)
    return 1/psq*(np.log(-psq/msq)+np.log(1-msq/psq))

# the numerical value of the taylor series to some order of the exact result
def fapprox(psq,msq,order):
    msq = complex(msq)
    s = 0
    for j in range(1,1+order):
        s -= pow(msq/psq,j)/j
    return (s+np.log(-psq/msq))/psq

if __name__ == "__main__":

    psq, msq = 4, 0.002
    name = "bubble1L_dotted_z"
    real_parameters = [psq,msq,1] # z set to 1

    # load c++ library
    intlib = IntegralLibrary(f"{name}/{name}_pylink.so".format(name))
    intlib.use_Qmc(transform="korobov3", fitfunction="polysingular")

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = intlib(real_parameters)

    # convert the result to sympy expressions
    result, error = map(sp.sympify, series_to_sympy(str_integral_with_prefactor))

    # access and print individual terms of the expansion
    print("Numerical Result")
    for power in [-2, -1, 0]:
        val = complex(result.coeff("eps", power))
        err = complex(error.coeff("eps", power))
        print(f"eps^{power:<2} {val: .5f} +/- {err:.5e}")

    order = 1
    print("psq, msq: {}, {}".format(psq,msq))
    print("exact result: {:.5f}".format(f(psq,msq)))
    print("approximate result to order {}: {:.5f}".format(order,fapprox(psq,msq,order)))
