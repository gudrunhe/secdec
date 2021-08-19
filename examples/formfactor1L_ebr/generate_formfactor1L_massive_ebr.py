#!/usr/bin/env python3

import pySecDec as psd

if __name__ == "__main__":

    # Example 6 of Bernd Jantzen (arXiv:1111.2589)
    
    # replace Qsq by one and add the Qsq prefactor later to factorize Qsq out
    li = psd.LoopIntegralFromGraph(
        internal_lines = [['0',[1,3]],['0',[2,3]],['m',[1,2]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        replacement_rules = [
                            ('p1*p1', '0'),
                            ('p2*p2', '0'),
                            ('p3*p3', '-Qsq'),
                            ('m**2', 't*Qsq'),
                            ('Qsq', '1'),
                            ]
    )
                        
    # find the regions
    generators_args = psd.loop_regions(
        name = 'formfactor1L_massive_ebr',
        loop_integral=li,
        smallness_parameter = 't',
        expansion_by_regions_order=0,
        decomposition_method = 'geometric',
        additional_prefactor=f"Qsq**({li.exponent_F})",
        add_monomial_regulator_power="n",
        )

    # write the code to sum up the regions
    psd.sum_package('formfactor1L_massive_ebr',
                generators_args,
                li.regulators,
                requested_orders = [0,0],
                real_parameters = ['Qsq','t'],
                complex_parameters = [])
