#!/usr/bin/env python3

import pySecDec as psd

if __name__ == "__main__":

    # Section 3.1 (Method of Regions for the Sudakov form factor) of Thomas Becher (arXiv:1803.04310)
    li = psd.LoopIntegralFromGraph(
        internal_lines = [['0',[1,3]],['0',[2,3]],['0',[1,2]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        replacement_rules = [
                            ('p3*p3', '-qsq'),
                            ('p1*p1', '-z*psq'),
                            ('p2*p2', '-z*lsq'),
                            ]
    )

    # find the regions
    generators_args = psd.loop_regions(
        name = 'formfactor1L_massless_ebr',
        loop_integral=li,
        smallness_parameter = 'z',
        expansion_by_regions_order=0)

    # write the code to sum up the regions
    psd.sum_package('formfactor1L_massless_ebr',
                generators_args,
                li.regulators,
                requested_orders = [0],
                real_parameters = ['qsq','lsq','psq','z'],
                complex_parameters = [])
