#!/usr/bin/env python3

import pySecDec as psd

if __name__ == "__main__":

    # Example 6 of Bernd Jantzen (arXiv:1111.2589)
    
    
    li = psd.LoopIntegralFromGraph(
        internal_lines = [['0',[1,3]],['0',[2,3]],['m',[1,2]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        powerlist=['1+n/2','1+n/3','1'],
        regulators=['eps','n'],
        replacement_rules = [
                            ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '-qsq'),
                            ('m**2', 'z*msq')
                            ]
    )
                        
    # find the regions
    generators_args = psd.loop_regions(
        name = 'formfactor1L_massive_ebr',
        loop_integral=li,
        smallness_parameter = 'z',
        expansion_by_regions_order=0,
        decomposition_method = 'geometric')

    # write the code to sum up the regions
    psd.sum_package('formfactor1L_massive_ebr',
                generators_args,
                li.regulators,
                requested_orders = [0,0],
                real_parameters = ['qsq','msq','z'],
                complex_parameters = [])
