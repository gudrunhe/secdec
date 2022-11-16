#!/usr/bin/env python3

from pySecDec import make_package

if __name__ == "__main__":

    make_package(
        name = 'easy',
        integration_variables = ['x','y'],
        regulators = ['eps'],

        requested_orders = [0],
        polynomials_to_decompose = ['(x+y)^(-2+eps)']
    )
