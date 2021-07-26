#!/usr/bin/env python3

# This example is the nonplanar two loop nonplanar example in section 3.3 in arXiv:2003.02451
# for example when using qsq, msq = -2.7, 0.2, all the poles cancel and only the finite term remains:
# (2.81872,-2.70076e-07) +/- (6.44647e-05,2.45386e-07)
from pySecDec.loop_integral.loop_regions import loop_regions
import pySecDec as psd
from pySecDec.algebra import Polynomial

import sympy as sp
import numpy as np

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['m',[1,2]],['m',[2,3]],['m',[3,4]],['m',[4,1]],[0,[1,5]],[0,[5,3]],],
    external_lines = [['p1',2],['p2',4],['q',5],],

    powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5","1+n1/7","1+n1/11"],
    # regulators=["eps"],
    regulators=["n1","eps"],
    Feynman_parameters=["x%i" % i for i in range(1,7)], # this renames the parameters

    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p1*p2', 'q**2/2'),
                            ('m**2', 'msq'),
                            ('q**2', 'qsq'),
                        ]
    )

    name = "nonplanar_2L"

    generators_args = psd.loop_integral.loop_regions(
        name = name,
        loop_integral=li,
        smallness_parameter = "msq",
        expansion_by_regions_order=0,
        form_work_space="5000M",
        enforce_complex=True,
        processes=None)#[4:5]


    for args in generators_args:
        # this removes the contour deformation and sign checking
        # useful to be able to run Euclidean points fast
        args.pop("contour_deformation_polynomial")
        args.pop("positive_polynomials")

    psd.code_writer.sum_package(name, [psd.make_package]*len(generators_args), generators_args, li.regulators,
                requested_orders = [0,0],
                real_parameters = ['qsq','msq'], complex_parameters = [])
