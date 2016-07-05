#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [['m',[3,4]],['m',[4,5]],['m',[3,5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]],
external_lines = [['p1',1],['p2',2],['p3',3]],

replacement_rules = [
                        ('p1*p1', 0),
                        ('p2*p2', 0),
                        ('p3*p3', 's'),
                        ('p4*p4', 0),
                        ('p1*p2', 's/2'),
                        ('p2*p3', '-s/2'),
                        ('p1*p3', '-s/2'),
                        ('m**2', 'msq')
                    ]
)


Mandelstam_symbols = ['s']
mass_symbols = ['msq']


loop_package(

target_directory = 'pySecDec_loop_integrals',

name = 'triangle',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '10G',

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
contour_deformation = True,

)
