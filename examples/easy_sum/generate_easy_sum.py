#!/usr/bin/env python3
from pySecDec import MakePackage
from pySecDec import sum_package

if __name__ == "__main__":

    coefficients = {
        "sum1": [
            '2*s',   # easy1
            '3*s'    # easy2
        ],
        "sum2": [
            's/(2*eps)', # easy1
            's*eps/3'  # easy2
        ]
    }

    integrals = [
        MakePackage('easy1',
            integration_variables = ['x','y'],
            polynomials_to_decompose = ['(x+y)^(-2+eps)'],
            ),
        MakePackage('easy2',
            integration_variables = ['x','y'],
            polynomials_to_decompose = ['(2*x+3*y)^(-1+eps)'],
            )
    ]
    
    # generate code sum of (int * coeff)
    sum_package('easy_sum', integrals,
        coefficients = coefficients, real_parameters=['s'],
        regulators=['eps'], requested_orders=[0])
