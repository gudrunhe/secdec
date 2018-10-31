#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd


li = psd.loop_integral.LoopIntegralFromGraph(

internal_lines = [ [0,[1,2]],[0,[1,3]],[0,[2,3]],[0,[2,5]],[0,[5,6]],[0,[6,3]],[0,[6,7]],[0,[5,7]],[0,[7,4]],[0,[4,7]],[0,[4,3]],[0,[4,1]] ],
external_lines = [['p',1],['p',2]],

replacement_rules = [ ('p*p','-1') ]

)

loop_package(

name = 'bubble6L',

loop_integral = li,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 4,

# the WorkSpace parameter for FORM
form_work_space = '2G',

# the method to be used for the sector decomposition
# valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
decomposition_method = 'geometric',
# if you choose ``geometric[_ku]`` and 'normaliz' is not in your
# $PATH, you can set the path to the 'normaliz' command-line
# executable here
#normaliz_executable='/path/to/normaliz',

# whether or not to produce code to perform the contour deformation
# contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
contour_deformation = False,

)

