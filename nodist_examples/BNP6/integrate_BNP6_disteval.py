from __future__ import print_function
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    integral = DistevalLibrary('BNP6/disteval/BNP6.json')

    # integrate
    result = integral(real_parameters=[9.,-2.5])

    print('Numerical Result')
    print(result)

    print('Analytic Result at (s,t,u)=(9,-2.5,-6.5)')
    print('eps^-2 with prefactor:', '0.10907856854318447 - 0.5799863360473464*I')
    print('eps^-1 with prefactor:', '-0.8876663743916553 + 4.360251717854891*I')
    print('eps^0 with prefactor: ', '0.7966721383373115 - 18.22048104236002*I')
