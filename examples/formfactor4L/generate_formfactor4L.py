#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import sympy as sp
import pySecDec as psd

li = psd.loop_integral.LoopIntegralFromGraph(
internal_lines = [[0,[1,4]],[0,[4,2]],[0,[2,5]],[0,[5,6]],[0,[6,7]],[0,[7,3]],[0,[3,8]],[0,[8,4]],[0,[1,9]],[0,[9,7]],[0,[9,5]],[0,[6,8]]],
external_lines = [['p1',1],['p2',2],['p3',3]],
powerlist=[2,1,1,1,1,1,1,1,1,2,1,1],
dimensionality='6-2*eps',
replacement_rules = [
                        ('p1*p1',-1),
                        ('p2*p2',0),
                        ('p3*p3','0'),
                        ('p1*p2','1/2'),
                        ('p2*p3','-1/2'),
                        ('p1*p3','1/2')
                    ]
)


Mandelstam_symbols = []
mass_symbols = []


loop_package(

name = 'formfactor4L',

loop_integral = li,

# normalization as defined in equation (2.4) of arXiv:1510.06758
additional_prefactor = sp.gamma(li.dimensionality/2 - 1) ** li.L,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 1,

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 4,

# the WorkSpace parameter for FORM
form_work_space = '1G',

# the method to be used for the sector decomposition
# valid values are ``iterative`` and ``geometric``
decomposition_method = 'geometric',
# if you choose ``geometric`` and 'normaliz' is not in your
# $PATH, you can set the path to the 'normaliz' command-line
# executable here
#normaliz_executable='/path/to/normaliz',

contour_deformation = False,

# there are no symmetries
use_dreadnaut = False,
use_Pak = False,

)

# analytic result available from arXiv:1510.06758, equation (7.1)
