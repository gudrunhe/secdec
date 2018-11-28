#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd


li = psd.loop_integral.LoopIntegralFromPropagators(

loop_momenta = ['k1','k2','k3'],
external_momenta = ['p'],

propagators = ['(k1+p)**2-m1sq','(k1-k3)**2-m1sq','k2**2-m2sq','(k2-k3)**2-m3sq'],
powerlist = [1,1,1,1],

replacement_rules = [ ('p*p','s') ]

)


Mandelstam_symbols = ['s']
mass_symbols = ['m1sq','m2sq','m3sq']


loop_package(

name = 'banana_3mass',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3)
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
contour_deformation = True,

)

