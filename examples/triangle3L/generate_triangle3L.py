#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [ [0,[1,2]], [0,[1,4]], [0,[1,5]], [0,[2,4]], [0,[2,5]], [0,[3,4]], [0,[3,5]] ],
external_lines = [['p1',1],['p2',2],['p3',3]],
powerlist=[1,1,1,1,1,1,1],

replacement_rules = [
                        ('p1*p1',0),
                        ('p2*p2',0),
                        ('p3*p3','s'),
                        ('p1*p2','s/2'),
                        ('p2*p3','-s/2'),
                        ('p1*p3','-s/2'),
			('s',-1)
                    ]
)


Mandelstam_symbols = []
mass_symbols = []


loop_package(

name = 'triangle3L',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

#additional_prefactor = 'gamma( (4-2*eps)/2-1 )**3',
additional_prefactor = '(-eps*gamma(-eps))**3', # work around slow gamma expansion in sympy 1.0

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 4,

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '100M',

# the method to be used for the sector decomposition
# valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
decomposition_method = 'iterative',
# if you choose ``geometric[_ku]`` and 'normaliz' is not in your
# $PATH, you can set the path to the 'normaliz' command-line
# executable here
#normaliz_executable='/path/to/normaliz',

contour_deformation = False

)
