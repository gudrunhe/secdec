#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [['m',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

replacement_rules = [
                        ('p1*p1', 's1'),
                        ('p2*p2', 0),
                        ('p3*p3', 0),
                        ('p4*p4', 0),
                        ('p3*p2', 't/2'),
                        ('p1*p2', 's/2-s1/2'),
                        ('p1*p4', 't/2-s1/2'),
                        ('p2*p4', 's1/2-t/2-s/2'),
                        ('p3*p4', 's/2'),
                        ('m**2', 'msq')
                    ]
)

Mandelstam_symbols = ['s','t','s1']
mass_symbols = ['msq']

loop_package(

name = 'box1L',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,

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

)
