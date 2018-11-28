#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd


li = psd.loop_integral.LoopIntegralFromGraph(

    internal_lines = [[0,[1,2]],[0,[2,6]],[0,[6,3]],[0,[3,4]],[0,[4,5]],[0,[5,7]],[0,[7,1]],[0,[7,6]]],
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4],['p5',5]],
    dimensionality='6-2*eps',
    replacement_rules = [
                            ('p1*p1', 0),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
                            ('p5*p5', 0),
                            ('p1*p2', 's12/2'),
                            ('p1*p3', '(s45-s12-s23)/2'),
                            ('p1*p4', '(s23-s51-s45)/2'),
                            ('p1*p5', 's51/2'),
                            ('p2*p3', 's23/2'),
                            ('p2*p4', '(-s23-s34+s51)/2'),
                            ('p2*p5', '(s34-s12-s51)/2'),
                            ('p3*p4', 's34/2'),
                            ('p3*p5', '(s12-s34-s45)/2'),
                            ('p4*p5', 's45/2'),
                        ]
)


Mandelstam_symbols = ['s12', 's23','s34','s45', 's51']
mass_symbols = []


loop_package(

name = 'pentabox_fin',

loop_integral = li,

real_parameters = Mandelstam_symbols + mass_symbols,

# the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
requested_order = 0,
additional_prefactor = '1/gamma(2*(2 + eps))',
# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 4,

# the WorkSpace parameter for FORM
form_work_space = '1G',

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

# no symmetries --> no need to run the full symmetry finder
use_Pak = False,

)
