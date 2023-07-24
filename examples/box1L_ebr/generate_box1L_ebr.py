#!/usr/bin/env python3

from pySecDec import sum_package, loop_regions
import pySecDec as psd

# This example is the one-loop box example in Go Mishima's paper arXiv:1812.04373

if __name__ == "__main__":

    # here we define the Feynman diagram
    li = psd.loop_integral.LoopIntegralFromGraph(
    internal_lines = [['mt',[3,1]],['mt',[1,2]],['mt',[2,4]],['mt',[4,3]]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
    powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5"],
    regulators=["n1","eps"],
    Feynman_parameters=["x%i" % i for i in range(1,5)], # renames the parameters to get the same polynomials as in 1812.04373

    replacement_rules = [
                            # note that in those relations all momenta are incoming
                            # general relations:
                            ('p4', '-p1-p2-p3'),
                            ('p1*p1', 'm1sq'),
                            ('p2*p2', 'm2sq'),
                            ('p3*p3', 'm3sq'),
                            ('p1*p2', 's/2-(m1sq+m2sq)/2'),
                            ('p1*p3', 't/2-(m1sq+m3sq)/2'),
                            ('p2*p3', 'u/2-(m2sq+m3sq)/2'),
                            ('u', '(m1sq+m2sq+m3sq+m4sq)-s-t'),
                            # relations for our specific case:
                            ('mt**2', 'mtsq'),
                            ('m1sq',0),
                            ('m2sq',0),
                            ('m3sq','mHsq'),
                            ('m4sq','mHsq'),
                            ('mHsq', 0),
                        ])

    # find the regions
    generators_args = loop_regions(
        name = "box1L_ebr",
        loop_integral=li,
        smallness_parameter = "mtsq",
        expansion_by_regions_order=0)

    # write the code to sum up the regions
    sum_package("box1L_ebr",
                generators_args,
                li.regulators,
                requested_orders = [0,0],
                real_parameters = ['s','t','mtsq'],
                complex_parameters = [])
