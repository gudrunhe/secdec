#!/usr/bin/env python3

import pySecDec as psd

if __name__ == "__main__":

    # Section 3.1 (Method of Regions for the Sudakov form factor) of Thomas Becher (arXiv:1803.04310)
    # use relative dimensionless parameters Psqr and Lsqr, which should be roughly of order one
    # replace Qsq with one and add the Qsq prefactor later to factorize Qsq out
    li = psd.LoopIntegralFromGraph(
        internal_lines = [['0',[1,3]],['0',[2,3]],['0',[1,2]]],
        external_lines = [['p1',1],['p2',2],['p3',3]],
        replacement_rules = [
                            ('p3*p3', '-Qsq'),
                            ('p1*p1', '-t*Psqr*Qsq'),
                            ('p2*p2', '-t*Lsqr*Qsq'),
                            ("Qsq", 1),
                            ],
    )

    # find the regions
    generators_args = psd.loop_regions(
        name = 'formfactor1L_massless_ebr',
        loop_integral=li,
        smallness_parameter = 't',
        expansion_by_regions_order=0,
        additional_prefactor=f"Qsq**({li.exponent_F})",
        )

    # write the code to sum up the regions
    psd.sum_package('formfactor1L_massless_ebr',
                generators_args,
                li.regulators,
                requested_orders = [0],
                real_parameters = ['Qsq','t','Psqr','Lsqr'],
                complex_parameters = [])
