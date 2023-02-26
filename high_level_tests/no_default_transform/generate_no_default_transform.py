#!/usr/bin/env python3
import shutil
from pySecDec import make_package

if __name__ == "__main__":

    name = 'no_default_transform'

    make_package(

    name=name,
    integration_variables = ['x'],

    # the order here defines the order of the expansion
    regulators = ['eps'],

    # the highest orders of the final regulator expansion
    # the order here matches the order of ``regulators``
    requested_orders = [1],

    polynomials_to_decompose = ['(x)**(-1+eps)'],
    polynomial_names = ['p'],
    other_polynomials = ['p**-2'],
    pylink_qmc_transforms=['none', 'korobov2x1', 'sidi3']
    )
