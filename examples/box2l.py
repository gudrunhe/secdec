#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd


li = psd.loop_integral.LoopIntegralFromPropagators(
propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2','(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2','(k1+p3)**2'],
loop_momenta = ['k1','k2'],
powerlist = [1,1,1,1,1,1,1,0],
external_momenta = ['p1','p2','p3','p4'],

replacement_rules = [
                        ('p1*p1', 0),
                        ('p2*p2', 0),
                        ('p3*p3', 0),
                        ('p4*p4', 0),
                        ('p1*p2', 's/2'),
                        ('p2*p3', 't/2'),
                        ('p1*p3', '-s/2-t/2')
                    ]

)


Mandelstam_symbols = ['s', 't']
mass_symbols = []


loop_package(

name = 'box2l',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '100M',

# whether or not to increase numerical stability
# Note: This is very extensive concerning both - the algebra and the numerics.
#       It should only be set to ``True`` if numerical instabilities occur.
stabilize = False,

# the method to be used for the sector decomposition
# valid values are ``iterative`` and ``geometric``
decomposition_method = 'geometric',

# whether or not to produce code to perform the contour deformation
# if ``True``, it can still be deactivated in the "config.hpp"
# if ``False``, no code for the contour deformation is generated
contour_deformation = False,

)
