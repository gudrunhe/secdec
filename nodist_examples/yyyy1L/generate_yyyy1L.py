#! /usr/bin/env python
from pySecDec.loop_integral import LoopPackage
from pySecDec.code_writer import sum_package

# import integrals and coefficients
import integrals
import coefficients

# 4-photon amplitude M++--
# Amp4y_ppmm = -8*( 1 + (t**2+u**2)/s *Box6dim(t,u) +  (t-u)/s*(  BubbleD(u)-BubbleD(t) ) )

if __name__ == '__main__':

    # define input to loop_package
    package_generators = []
    names = ['bubble_u', 'bubble_t', 'box_6', 'box_8']
    for i in range(len(integrals.I)):
        package_generators.append(
            LoopPackage(
                name = names[i],
                loop_integral = integrals.I[i],
                requested_orders = [0]
                )
        )

    # generate code that will calculate the sum of all integrals times the corresponding coefficients
    sum_package(
        'yyyy1L',
        package_generators,
        regulators = ['eps'],
        requested_orders = [0],
        real_parameters = ['t','u'],
        coefficients = coefficients.coeff,
        complex_parameters = []
    )

