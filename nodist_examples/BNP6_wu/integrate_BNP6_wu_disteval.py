from __future__ import print_function
from pySecDec.integral_interface import DistevalLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    bnp6_wu = DistevalLibrary('BNP6_wu/disteval/BNP6_wu.json')

    # integrate
    result = bnp6_wu(parameters={"s": 9., "t": -2.5, "a0": 30./40., "a1": 42./107., "a2": 51./65., "a3": 67./89., "a4": 79./55., "a5": 88./33.})

    print('Numerical Result')
    print(result)

    print('Analytic Result at (s,t,u)=(9,-2.5,-6.5)')
    print('eps^-2 with prefactor:', '0.10907856854318447 - 0.5799863360473464*I')
    print('eps^-1 with prefactor:', '-0.8876663743916553 + 4.360251717854891*I')
    print('eps^0 with prefactor: ', '0.7966721383373115 - 18.22048104236002*I')
