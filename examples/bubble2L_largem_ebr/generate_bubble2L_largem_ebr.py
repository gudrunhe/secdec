#!/usr/bin/env python3

from pySecDec.loop_integral import loop_regions
from pySecDec.code_writer import sum_package, make_package
import pySecDec as psd

import pySecDec as psd

# example is I3tilde of hep-ph/9605392

if __name__ == "__main__":

    # define Feynman Integral
    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['k1**2-msq_','(k1+p1)**2-msq_','(k1-k2)**2-msq_','k2**2','(k2+p1)**2'],
    loop_momenta = ['k1','k2'],
    external_momenta = ['p1'],
    replacement_rules = [('p1*p1', 'psq'), ('msq_','msq/z')])

    # find the regions and expand the integrals using expansion by regions
    regions_generator_args = loop_regions(
        name = 'bubble2L_largem_ebr',
        loop_integral = li,
        smallness_parameter = 'z',
        additional_prefactor = '-1',
        expansion_by_regions_order = 1)

    # generate code that will calculate the sum of all regions and all orders in
    # the smallness parameter
    sum_package('bubble2L_largem_ebr',
        regions_generator_args, regulators = ['eps'],requested_orders = [0],
        real_parameters = ['psq','msq','z'],
        complex_parameters = [])
