#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromPropagators(

loop_momenta = ['k1','k2'],
external_momenta = ['p1','p2','p3'],

propagators = ['(k1+p1)**2', '(k2+p1+p2)**2', '(k2-p3)**2', '(k1-k2+p3)**2', '(k1-k2)**2-mm'],
powerlist = [1, 1, 1, 1, 1],

replacement_rules = [
                        ('p1*p1', '0'),
                        ('p2*p2', '0'),
                        ('p3*p3', '0'),
                        ('p1*p2', 's/2'),
                        ('p1*p3', 't/2'),
                        ('p2*p3', '(mm-s-t)/2')
                    ]

)

Mandelstam_symbols = ['s','t']
mass_symbols = ['mm']

loop_package(

name = 'Nbox2L_split_c',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

additional_prefactor = '1',

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 4,

# the WorkSpace parameter for FORM
form_work_space = '1G',

# the method to be used for the sector decomposition
# valid values are ``iterative`` and ``geometric``
decomposition_method = 'iterative',
             
contour_deformation=True,

# there are singularities at one due to ``p4*p4 = mm``
#split = True # not needed for (quasi-)finite integrals

)

