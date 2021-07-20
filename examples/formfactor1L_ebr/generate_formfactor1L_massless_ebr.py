#!/usr/bin/env python3

from pySecDec import sum_package
from pySecDec.loop_integral import loop_regions
import pySecDec as psd

if __name__ == "__main__":

    # Section 3.1 (Method of Regions for the Sudakov form factor) of Thomas Becher (arXiv:1803.04310)
    li = psd.loop_integral.LoopIntegralFromGraph(
        internal_lines = [['0',[1,3]],['0',[2,3]],['0',[1,2]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        replacement_rules = [
                            ('p3*p3', '-qsq'),
                            ('p1*p1', '-z**2*psq'),
                            ('p2*p2', '-z**2*lsq'),
                            ]
    )
    
    # find the regions
    generators_args = loop_regions(
        name = 'formfactor1L_massless_ebr',
        loop_integral=li,
        smallness_parameter = 'z',
        expansion_by_regions_order=0)

    # write the code to sum up the regions
    sum_package('formfactor1L_massless_ebr',
                generators_args,
                li.regulators,
                requested_orders = [0],
                real_parameters = ['qsq','lsq','psq','z'],
                complex_parameters = [])
