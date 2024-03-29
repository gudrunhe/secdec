#!/usr/bin/env python3
import pySecDec as psd

if __name__ == "__main__":

    li = psd.LoopIntegralFromGraph(

    internal_lines = [ [0,[1,2]],[0,[1,3]],[0,[2,3]],[0,[2,5]],[0,[5,6]],[0,[6,3]],[0,[6,7]],[0,[5,7]],[0,[7,4]],[0,[4,7]],[0,[4,3]],[0,[4,1]] ],
    external_lines = [['p',1],['p',2]],

    replacement_rules = [ ('p*p','-1') ]

    )

    psd.loop_package(

    name = 'bubble6L',

    loop_integral = li,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_orders = [0],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 4,

    # the WorkSpace parameter for FORM
    form_work_space = '2G',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'geometric',

    # whether or not to produce code to perform the contour deformation
    # contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    contour_deformation = False,

    )

