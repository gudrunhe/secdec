#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromPropagators(
propagators = ['k1**2-msq','(k1+p1+p2)**2-msq','k2**2-msq','(k2+p1+p2)**2-msq','(k1+p1)**2-msq','(k1-k2)**2','(k2-p3)**2-msq','(k2+p1)**2','(k1-p3)**2'],
powerlist = [1,1,0,1,1,1,1,0,0],
loop_momenta = ['k1','k2'],
replacement_rules = [
                        ('p1*p1',0),
                        ('p2*p2',0),
                        ('p3*p3',0),
                        ('p1*p2','s/2'),
                        ('p2*p3','pp4/2-s/2-t/2'),
                        ('p1*p3','t/2')
                    ]
)


Mandelstam_symbols = ['s','t','pp4']
mass_symbols = ['msq']


loop_package(

name = 'elliptic2L_physical',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

#additional_prefactor = '(s/msq)**(3/2)',

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '100M',

# the method to be used for the sector decomposition
# valid values are ``iterative`` and ``geometric``
decomposition_method = 'iterative',
# if you choose ``geometric`` and 'normaliz' is not in your
# $PATH, you can set the path to the 'normaliz' command-line
# executable here
#normaliz_executable='/path/to/normaliz',

contour_deformation = True

)
