#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromPropagators(
loop_momenta = ['k1'],
external_momenta = ['p1','p2','p3','p4'],
Lorentz_indices = ['mu1','mu2','mu3','mu4'],

propagators = ['k1**2-msq','(k1+p1)**2-msq','(k1+p1+p2)**2-msq','(k1+p1+p2+p3)**2-msq'],
powerlist = [1,1,1,1],
numerator = 'k1(mu1)*p1(mu1)*k1(mu2)*p2(mu2)*k1(mu3)*p3(mu3)*k1(mu4)*p4(mu4)',

replacement_rules = [
                        ('p1*p1', 0),
                        ('p2*p2', 0),
                        ('p3*p3', 0),
                        ('p4*p4', 0),
                        ('p3*p1', '-t/2-s/2'),
                        ('p3*p2', 't/2'),
                        ('p1*p2', 's/2'),
                        ('p1*p4', 't/2'),
                        ('p2*p4', '-t/2-s/2'),
                        ('p3*p4', 's/2'),
                        ('m**2', 'msq')
                    ]
)

Mandelstam_symbols = ['s','t']
mass_symbols = ['msq']

loop_package(

name = 'box1L_rank4',
#name = 'test',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

# the optimization level to use in FORM (can be 0, 1, 2, 3)
form_optimization_level = 3,

# the WorkSpace parameter for FORM
form_work_space = '100M',

# the method to be used for the sector decomposition
# valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
decomposition_method = 'iterative',
# if you choose ``geometric[_ku]`` and 'normaliz' is not in your
# $PATH, you can set the path to the 'normaliz' command-line
# executable here
#normaliz_executable='/path/to/normaliz',
contour_deformation = True,

)
