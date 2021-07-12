#! /usr/bin/env python3
from pySecDec.loop_integral import LoopPackage
from pySecDec.code_writer import sum_package

# import integrals and coefficients
import integrals
import coefficients

# 4-photon amplitude M++--

if __name__ == '__main__':

    # define input to loop_package
    package_generators = []
    names = ['bubble_u', 'bubble_t', 'box_6', 'box_8']
    for name, integral in zip(names, integrals.I):
        package_generators.append(
            LoopPackage(
                name = name,
                loop_integral = integral,
                requested_orders = [0]
                )
        )

    # generate code sum of (int * coeff)
    sum_package(
        'yyyy1L',
        package_generators,
        regulators = ['eps'],
        requested_orders = [0],
        real_parameters = ['t','u'],
        coefficients = coefficients.coeff,
        complex_parameters = []
    )

