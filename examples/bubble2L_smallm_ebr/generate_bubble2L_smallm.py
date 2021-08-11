#!/usr/bin/env python3

import pySecDec as psd

# example is 2-loop bubble with five propagators, two of them massive

if __name__ == "__main__":

    # define Feynman Integral
    li = psd.LoopIntegralFromPropagators(
    propagators = ['k1**2-msq_','(k1+k2)**2-msq_','(k1+p1)**2','k2**2','(k1+k2+p1)**2'],
    loop_momenta = ['k1','k2'],
    external_momenta = ['p1'],
    replacement_rules = [('p1*p1', 'psq'), ('msq_','z*msq')])

    # find the regions and expand the integrals using expansion by regions
    regions_generator_args = psd.loop_regions(
        name = 'bubble2L_smallm',
        loop_integral = li,
        smallness_parameter = 'z',
        expansion_by_regions_order = 2) # this has to be 2 to catch the ultra-soft region

    # generate code that will calculate the sum of all regions and all orders in
    # the smallness parameter
    psd.sum_package('bubble2L_smallm',
        regions_generator_args, regulators = ['eps'],requested_orders = [0],
        real_parameters = ['psq','msq','z'],
        complex_parameters = [])
