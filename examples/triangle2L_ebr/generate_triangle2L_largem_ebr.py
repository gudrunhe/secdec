#!/usr/bin/env python3

from pySecDec import LoopIntegralFromPropagators
from pySecDec import loop_regions
from pySecDec import sum_package

if __name__ == "__main__":

    # this object represents the Feynman graph
    li = LoopIntegralFromPropagators(
            propagators = ['(k1+p1)**2','(k1-p2)**2','(k2+p1)**2','(k2-p2)**2','k2**2','(k1-k2)**2-msq_'],
            loop_momenta = ['k1','k2'],
            external_momenta = ['p1','p2'],
                replacement_rules = [
                                        ('p1*p1', 0),
                                        ('p2*p2', 0),
                                        ('p1*p2', 's/2'),
                                        ('msq_', 'msq/z')
                                    ]
    )
    
    # find the regions
    generators_args = loop_regions(
        name = 'triangle2L_largem_ebr',
        loop_integral = li,
        smallness_parameter = 'z',
        expansion_by_regions_order = 1,
        additional_prefactor = '1/gamma(2+2*eps)',
        decomposition_method = 'geometric')

    # write the code to sum up the regions
    sum_package('triangle2L_largem_ebr',
                generators_args,
                li.regulators,
                requested_orders = [0],
                real_parameters = ['s','msq','z'],
                complex_parameters = [])
