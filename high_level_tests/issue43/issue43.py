#!/usr/bin/env python3
import shutil
from pySecDec import make_package

if __name__ == "__main__":

    name = 'issue43'

    make_package(

    name=name,
    integration_variables = ['z'],
    prefactor = '1',
    regulators = ['eps'],
    requested_orders = [0],

    decomposition_method = 'geometric_no_primary', # or iterative_no_primary
    polynomials_to_decompose = ['(1-z)^(-1-eps)'],
    split=True,

    )
