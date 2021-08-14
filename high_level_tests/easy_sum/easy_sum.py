#!/usr/bin/env python3
from pySecDec.code_writer.sum_package import Coefficient
from pySecDec.make_package import MakePackage
from pySecDec.code_writer import sum_package

if __name__ == "__main__":

    common_args = {}
    common_args['real_parameters'] = ['s']
    common_args['regulators'] = ['eps']
    common_args['requested_orders'] = [0]

    coefficients = [
        [ # sum1
            Coefficient(['2*s'],['1'],['s']),   # easy1
            Coefficient(['3*s'],['1'],['s'])    # easy2
        ],
        [ # sum2
            Coefficient(['s'],['2*eps'],['s']), # easy1
            Coefficient(['s*eps'],['3'],['s'])  # easy2
        ]
    ]

    integrals = [
        MakePackage('easy1',
            integration_variables = ['x','y'],
            polynomials_to_decompose = ['(x+y)^(-2+eps)'],
            **common_args),
        MakePackage('easy2',
            integration_variables = ['x','y'],
            polynomials_to_decompose = ['(2*x+3*y)^(-1+eps)'],
            **common_args)
    ]
    
    # generate code sum of (int * coeff)
    sum_package('easy_sum', integrals,
        coefficients = coefficients, **common_args)
