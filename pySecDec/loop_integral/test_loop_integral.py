from . import *
from ..misc import sympify_expression
import numpy as np
import sympy as sp
import os
import shutil
import tempfile
import unittest
from math import floor
from nose.plugins.attrib import attr

def uf_from_propagators_generic(test_case, loop_momenta, propagators, result_u, result_f):
    loop_integral = LoopIntegralFromPropagators(propagators, loop_momenta, Feynman_parameters='Feynman')

    u,f = loop_integral.U, loop_integral.F

    sympy_u = sympify_expression(str(u))
    sympy_f = sympify_expression(str(f))

    result_u = sympify_expression(str(result_u).replace('x','Feynman'))
    result_f = sympify_expression(str(result_f).replace('x','Feynman'))

    zerou = (sympy_u - result_u).expand()
    zerof = (sympy_f - result_f).expand()

    test_case.assertEqual(zerou,0)
    test_case.assertEqual(zerof,0)

    expo_u, expo_f = loop_integral.exponent_U, loop_integral.exponent_F
    test_case.assertEqual(expo_u, len(propagators) - sympify_expression('2-eps') * (1 + len(loop_momenta)))
    test_case.assertEqual(expo_f, -(len(propagators) - sympify_expression('2-eps') * len(loop_momenta)))

#@attr('active')
class TestUF_FromPropagators(unittest.TestCase):
    def test_massive(self):
        # from SecDec: loop -> demos -> 6_geomethod_2L
        uf_from_propagators_generic(self,
            loop_momenta = ['k','l'],
            propagators = ['k**2','(k-p1)**2-m**2','(k+p2)**2-m**2','(k-l)**2',
                           '(l-p1)**2-m**2','(l+p2)**2-m**2','(l+p2+p3)**2'],
            result_u = """x3*(x4 + x5 + x6) + x0*(x3 + x4 + x5 + x6) +
                          x1*(x3 + x4 + x5 + x6) + x2*(x3 + x4 + x5 + x6)""",
            result_f = """m**2*x0*x1*x3 + m**2*x0*x1*x4 + m**2*x0*x1*x5 + m**2*x0*x1*x6 + m**2*x0*x2*x3 +
                          m**2*x0*x2*x4 + m**2*x0*x2*x5 + m**2*x0*x2*x6 + m**2*x0*x3*x4 + m**2*x0*x3*x5 +
                          m**2*x0*x4**2 + 2*m**2*x0*x4*x5 + m**2*x0*x4*x6 + m**2*x0*x5**2 + m**2*x0*x5*x6 +
                          m**2*x1**2*x3 + m**2*x1**2*x4 + m**2*x1**2*x5 + m**2*x1**2*x6 + 2*m**2*x1*x2*x3 +
                          2*m**2*x1*x2*x4 + 2*m**2*x1*x2*x5 + 2*m**2*x1*x2*x6 + 2*m**2*x1*x3*x4 +
                          2*m**2*x1*x3*x5 + m**2*x1*x3*x6 + m**2*x1*x4**2 + 2*m**2*x1*x4*x5 + m**2*x1*x4*x6 +
                          m**2*x1*x5**2 + m**2*x1*x5*x6 + m**2*x2**2*x3 + m**2*x2**2*x4 + m**2*x2**2*x5 +
                          m**2*x2**2*x6 + 2*m**2*x2*x3*x4 + 2*m**2*x2*x3*x5 + m**2*x2*x3*x6 + m**2*x2*x4**2 +
                          2*m**2*x2*x4*x5 + m**2*x2*x4*x6 + m**2*x2*x5**2 + m**2*x2*x5*x6 + m**2*x3*x4**2 +
                          2*m**2*x3*x4*x5 + m**2*x3*x4*x6 + m**2*x3*x5**2 + m**2*x3*x5*x6 - p1**2*x0*x1*x3 -
                          p1**2*x0*x1*x4 - p1**2*x0*x1*x5 - p1**2*x0*x1*x6 - p1**2*x0*x3*x4 - p1**2*x0*x4*x5 -
                          p1**2*x0*x4*x6 - p1**2*x1*x2*x3 - p1**2*x1*x2*x4 - p1**2*x1*x2*x5 - p1**2*x1*x2*x6 -
                          p1**2*x1*x3*x5 - p1**2*x1*x3*x6 - p1**2*x1*x4*x5 - p1**2*x1*x4*x6 - p1**2*x2*x3*x4 -
                          p1**2*x2*x4*x5 - p1**2*x2*x4*x6 - p1**2*x3*x4*x5 - p1**2*x3*x4*x6 - 2*p1*p2*x0*x4*x5 -
                          2*p1*p2*x0*x4*x6 - 2*p1*p2*x1*x2*x3 - 2*p1*p2*x1*x2*x4 - 2*p1*p2*x1*x2*x5 -
                          2*p1*p2*x1*x2*x6 - 2*p1*p2*x1*x3*x5 - 2*p1*p2*x1*x3*x6 - 2*p1*p2*x1*x4*x5 -
                          2*p1*p2*x1*x4*x6 - 2*p1*p2*x2*x3*x4 - 2*p1*p2*x2*x4*x5 - 2*p1*p2*x2*x4*x6 -
                          2*p1*p2*x3*x4*x5 - 2*p1*p2*x3*x4*x6 - 2*p1*p3*x0*x4*x6 - 2*p1*p3*x1*x3*x6 -
                          2*p1*p3*x1*x4*x6 - 2*p1*p3*x2*x4*x6 - 2*p1*p3*x3*x4*x6 - p2**2*x0*x2*x3 -
                          p2**2*x0*x2*x4 - p2**2*x0*x2*x5 - p2**2*x0*x2*x6 - p2**2*x0*x3*x5 - p2**2*x0*x3*x6 -
                          p2**2*x0*x4*x5 - p2**2*x0*x4*x6 - p2**2*x1*x2*x3 - p2**2*x1*x2*x4 - p2**2*x1*x2*x5 -
                          p2**2*x1*x2*x6 - p2**2*x1*x3*x5 - p2**2*x1*x3*x6 - p2**2*x1*x4*x5 - p2**2*x1*x4*x6 -
                          p2**2*x2*x3*x4 - p2**2*x2*x4*x5 - p2**2*x2*x4*x6 - p2**2*x3*x4*x5 - p2**2*x3*x4*x6 -
                          2*p2*p3*x0*x3*x6 - 2*p2*p3*x0*x4*x6 - 2*p2*p3*x1*x3*x6 - 2*p2*p3*x1*x4*x6 -
                          2*p2*p3*x2*x4*x6 - 2*p2*p3*x3*x4*x6 - p3**2*x0*x3*x6 - p3**2*x0*x4*x6 -
                          p3**2*x0*x5*x6 - p3**2*x1*x3*x6 - p3**2*x1*x4*x6 - p3**2*x1*x5*x6 -
                          p3**2*x2*x3*x6 - p3**2*x2*x4*x6 - p3**2*x2*x5*x6 - p3**2*x3*x4*x6 -
                          p3**2*x3*x5*x6""")

    def test_bubble_1l(self):
        uf_from_propagators_generic(self,
            loop_momenta = ['k1'],
            propagators = ['k1**2','(k1-p1)**2'],
            result_u = "x0 + x1",
            result_f = "-p1**2*x0*x1")

    def test_box_2l(self):
        uf_from_propagators_generic(self,
            loop_momenta = ['k1','k2'],
            propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2',
                           '(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2'],
            result_u = """x0*x3 + x0*x4 + x0*x5 + x0*x6 + x1*x3 +
                          x1*x4 + x1*x5 + x1*x6 + x2*x3 + x2*x4 +
                          x2*x5 + x2*x6 + x3*x4 + x3*x5 + x3*x6""",
            result_f = """-p1**2*x0*x2*x3 - p1**2*x0*x2*x4 - p1**2*x0*x2*x5 -
                           p1**2*x0*x2*x6 - p1**2*x0*x3*x5 - p1**2*x0*x4*x5 -
                           p1**2*x0*x5*x6 - p1**2*x1*x2*x3 - p1**2*x1*x2*x4 -
                           p1**2*x1*x2*x5 - p1**2*x1*x2*x6 - p1**2*x1*x3*x5 -
                           p1**2*x1*x4*x5 - p1**2*x1*x5*x6 - p1**2*x2*x3*x4 -
                           p1**2*x2*x3*x6 - p1**2*x2*x4*x5 - p1**2*x2*x5*x6 -
                           p1**2*x3*x4*x5 - p1**2*x3*x5*x6 - 2*p1*p2*x0*x4*x5 -
                           2*p1*p2*x0*x5*x6 - 2*p1*p2*x1*x2*x3 -
                           2*p1*p2*x1*x2*x4 - 2*p1*p2*x1*x2*x5 -
                           2*p1*p2*x1*x2*x6 - 2*p1*p2*x1*x3*x5 -
                           2*p1*p2*x1*x4*x5 - 2*p1*p2*x1*x5*x6 -
                           2*p1*p2*x2*x3*x4 - 2*p1*p2*x2*x3*x6 -
                           2*p1*p2*x2*x4*x5 - 2*p1*p2*x2*x5*x6 -
                           2*p1*p2*x3*x4*x5 - 2*p1*p2*x3*x5*x6 -
                           2*p1*p3*x0*x5*x6 - 2*p1*p3*x1*x5*x6 -
                           2*p1*p3*x2*x3*x6 - 2*p1*p3*x2*x5*x6 -
                           2*p1*p3*x3*x5*x6 - p2**2*x0*x1*x3 - p2**2*x0*x1*x4 -
                           p2**2*x0*x1*x5 - p2**2*x0*x1*x6 - p2**2*x0*x3*x4 -
                           p2**2*x0*x3*x6 - p2**2*x0*x4*x5 - p2**2*x0*x5*x6 -
                           p2**2*x1*x2*x3 - p2**2*x1*x2*x4 - p2**2*x1*x2*x5 -
                           p2**2*x1*x2*x6 - p2**2*x1*x3*x5 - p2**2*x1*x4*x5 -
                           p2**2*x1*x5*x6 - p2**2*x2*x3*x4 - p2**2*x2*x3*x6 -
                           p2**2*x2*x4*x5 - p2**2*x2*x5*x6 - p2**2*x3*x4*x5 -
                           p2**2*x3*x5*x6 - 2*p2*p3*x0*x3*x6 -
                           2*p2*p3*x0*x5*x6 - 2*p2*p3*x1*x5*x6 -
                           2*p2*p3*x2*x3*x6 - 2*p2*p3*x2*x5*x6 -
                           2*p2*p3*x3*x5*x6 - p3**2*x0*x3*x6 - p3**2*x0*x4*x6 -
                           p3**2*x0*x5*x6 - p3**2*x1*x3*x6 - p3**2*x1*x4*x6 -
                           p3**2*x1*x5*x6 - p3**2*x2*x3*x6 - p3**2*x2*x4*x6 -
                           p3**2*x2*x5*x6 - p3**2*x3*x4*x6 - p3**2*x3*x5*x6""")

    @attr('slow')
    def test_formfactor_4l(self):
        uf_from_propagators_generic(self,
            loop_momenta = ['k1','k2','k3','k4'],
            propagators = ['k1**2','k2**2','k3**2','(k3-p1)**2','(k2-k3)**2','(k2-k3-k4)**2',
                           '(k1-k2)**2','(k1-p1-p2)**2','(k1-k2+k3+k4-p1-p2)**2',
                           '(-k1+k2-k4+p2)**2','(-k4+p2)**2','k4**2'],
            result_u = """x0*x1*x2*x5 + x0*x1*x3*x5 + x0*x1*x4*x5 + x0*x2*x4*x5 + x0*x3*x4*x5 + x0*x2*x5*x6 + x1*x2*x5*x6 + x0*x3*x5*x6 + x1*x3*x5*x6 + x0*x4*x5*x6 +
                          x1*x4*x5*x6 + x2*x4*x5*x6 + x3*x4*x5*x6 + x1*x2*x5*x7 + x1*x3*x5*x7 + x1*x4*x5*x7 + x2*x4*x5*x7 + x3*x4*x5*x7 + x2*x5*x6*x7 + x3*x5*x6*x7 +
                          x4*x5*x6*x7 + x0*x1*x2*x8 + x0*x1*x3*x8 + x0*x1*x4*x8 + x0*x2*x4*x8 + x0*x3*x4*x8 + x1*x2*x5*x8 + x1*x3*x5*x8 + x1*x4*x5*x8 + x2*x4*x5*x8 +
                          x3*x4*x5*x8 + x0*x2*x6*x8 + x1*x2*x6*x8 + x0*x3*x6*x8 + x1*x3*x6*x8 + x0*x4*x6*x8 + x1*x4*x6*x8 + x2*x4*x6*x8 + x3*x4*x6*x8 + x2*x5*x6*x8 +
                          x3*x5*x6*x8 + x4*x5*x6*x8 + x1*x2*x7*x8 + x1*x3*x7*x8 + x1*x4*x7*x8 + x2*x4*x7*x8 + x3*x4*x7*x8 + x2*x6*x7*x8 + x3*x6*x7*x8 + x4*x6*x7*x8 +
                          x0*x1*x2*x9 + x0*x1*x3*x9 + x0*x1*x4*x9 + x0*x2*x4*x9 + x0*x3*x4*x9 + x0*x1*x5*x9 + x1*x2*x5*x9 + x1*x3*x5*x9 + x0*x4*x5*x9 + x1*x4*x5*x9 +
                          x2*x4*x5*x9 + x3*x4*x5*x9 + x0*x2*x6*x9 + x1*x2*x6*x9 + x0*x3*x6*x9 + x1*x3*x6*x9 + x0*x4*x6*x9 + x1*x4*x6*x9 + x2*x4*x6*x9 + x3*x4*x6*x9 +
                          x0*x5*x6*x9 + x1*x5*x6*x9 + x2*x5*x6*x9 + x3*x5*x6*x9 + x1*x2*x7*x9 + x1*x3*x7*x9 + x1*x4*x7*x9 + x2*x4*x7*x9 + x3*x4*x7*x9 + x1*x5*x7*x9 +
                          x4*x5*x7*x9 + x2*x6*x7*x9 + x3*x6*x7*x9 + x4*x6*x7*x9 + x5*x6*x7*x9 + x0*x1*x8*x9 + x0*x4*x8*x9 + x1*x5*x8*x9 + x4*x5*x8*x9 + x0*x6*x8*x9 +
                          x1*x6*x8*x9 + x4*x6*x8*x9 + x5*x6*x8*x9 + x1*x7*x8*x9 + x4*x7*x8*x9 + x6*x7*x8*x9 + x0*x1*x2*x10 + x0*x1*x3*x10 + x0*x1*x4*x10 + x0*x2*x4*x10 +
                          x0*x3*x4*x10 + x0*x1*x5*x10 + x0*x2*x5*x10 + x0*x3*x5*x10 + x0*x2*x6*x10 + x1*x2*x6*x10 + x0*x3*x6*x10 + x1*x3*x6*x10 + x0*x4*x6*x10 + x1*x4*x6*x10 +
                          x2*x4*x6*x10 + x3*x4*x6*x10 + x0*x5*x6*x10 + x1*x5*x6*x10 + x2*x5*x6*x10 + x3*x5*x6*x10 + x1*x2*x7*x10 + x1*x3*x7*x10 + x1*x4*x7*x10 + x2*x4*x7*x10 +
                          x3*x4*x7*x10 + x1*x5*x7*x10 + x2*x5*x7*x10 + x3*x5*x7*x10 + x2*x6*x7*x10 + x3*x6*x7*x10 + x4*x6*x7*x10 + x5*x6*x7*x10 + x0*x1*x8*x10 + x0*x2*x8*x10 +
                          x1*x2*x8*x10 + x0*x3*x8*x10 + x1*x3*x8*x10 + x1*x4*x8*x10 + x2*x4*x8*x10 + x3*x4*x8*x10 + x1*x5*x8*x10 + x2*x5*x8*x10 + x3*x5*x8*x10 + x0*x6*x8*x10 +
                          x1*x6*x8*x10 + x4*x6*x8*x10 + x5*x6*x8*x10 + x1*x7*x8*x10 + x2*x7*x8*x10 + x3*x7*x8*x10 + x6*x7*x8*x10 + x0*x2*x9*x10 + x1*x2*x9*x10 + x0*x3*x9*x10 +
                          x1*x3*x9*x10 + x0*x4*x9*x10 + x1*x4*x9*x10 + x2*x4*x9*x10 + x3*x4*x9*x10 + x0*x5*x9*x10 + x1*x5*x9*x10 + x2*x5*x9*x10 + x3*x5*x9*x10 + x2*x7*x9*x10 +
                          x3*x7*x9*x10 + x4*x7*x9*x10 + x5*x7*x9*x10 + x0*x8*x9*x10 + x1*x8*x9*x10 + x4*x8*x9*x10 + x5*x8*x9*x10 + x7*x8*x9*x10 + x0*x1*x2*x11 + x0*x1*x3*x11 +
                          x0*x1*x4*x11 + x0*x2*x4*x11 + x0*x3*x4*x11 + x0*x1*x5*x11 + x0*x2*x5*x11 + x0*x3*x5*x11 + x0*x2*x6*x11 + x1*x2*x6*x11 + x0*x3*x6*x11 + x1*x3*x6*x11 +
                          x0*x4*x6*x11 + x1*x4*x6*x11 + x2*x4*x6*x11 + x3*x4*x6*x11 + x0*x5*x6*x11 + x1*x5*x6*x11 + x2*x5*x6*x11 + x3*x5*x6*x11 + x1*x2*x7*x11 + x1*x3*x7*x11 +
                          x1*x4*x7*x11 + x2*x4*x7*x11 + x3*x4*x7*x11 + x1*x5*x7*x11 + x2*x5*x7*x11 + x3*x5*x7*x11 + x2*x6*x7*x11 + x3*x6*x7*x11 + x4*x6*x7*x11 + x5*x6*x7*x11 +
                          x0*x1*x8*x11 + x0*x2*x8*x11 + x1*x2*x8*x11 + x0*x3*x8*x11 + x1*x3*x8*x11 + x1*x4*x8*x11 + x2*x4*x8*x11 + x3*x4*x8*x11 + x1*x5*x8*x11 + x2*x5*x8*x11 +
                          x3*x5*x8*x11 + x0*x6*x8*x11 + x1*x6*x8*x11 + x4*x6*x8*x11 + x5*x6*x8*x11 + x1*x7*x8*x11 + x2*x7*x8*x11 + x3*x7*x8*x11 + x6*x7*x8*x11 + x0*x2*x9*x11 +
                          x1*x2*x9*x11 + x0*x3*x9*x11 + x1*x3*x9*x11 + x0*x4*x9*x11 + x1*x4*x9*x11 + x2*x4*x9*x11 + x3*x4*x9*x11 + x0*x5*x9*x11 + x1*x5*x9*x11 + x2*x5*x9*x11 +
                          x3*x5*x9*x11 + x2*x7*x9*x11 + x3*x7*x9*x11 + x4*x7*x9*x11 + x5*x7*x9*x11 + x0*x8*x9*x11 + x1*x8*x9*x11 + x4*x8*x9*x11 + x5*x8*x9*x11 + x7*x8*x9*x11""",
            result_f =    """-p1**2*x0*x1*x2*x3*x5 - p1**2*x0*x1*x3*x4*x5 - p1**2*x0*x2*x3*x4*x5 - p1**2*x0*x2*x3*x5*x6 - p1**2*x1*x2*x3*x5*x6 - p1**2*x0*x3*x4*x5*x6 - p1**2*x1*x3*x4*x5*x6 -
                              p1**2*x2*x3*x4*x5*x6 - p1**2*x0*x1*x2*x5*x7 - 2*p1*p2*x0*x1*x2*x5*x7 - p2**2*x0*x1*x2*x5*x7 - p1**2*x0*x1*x3*x5*x7 - 2*p1*p2*x0*x1*x3*x5*x7 - p2**2*x0*x1*x3*x5*x7 -
                              p1**2*x1*x2*x3*x5*x7 - p1**2*x0*x1*x4*x5*x7 - 2*p1*p2*x0*x1*x4*x5*x7 - p2**2*x0*x1*x4*x5*x7 - p1**2*x0*x2*x4*x5*x7 - 2*p1*p2*x0*x2*x4*x5*x7 - p2**2*x0*x2*x4*x5*x7 -
                              p1**2*x0*x3*x4*x5*x7 - 2*p1*p2*x0*x3*x4*x5*x7 - p2**2*x0*x3*x4*x5*x7 - p1**2*x1*x3*x4*x5*x7 - p1**2*x2*x3*x4*x5*x7 - p1**2*x0*x2*x5*x6*x7 - 2*p1*p2*x0*x2*x5*x6*x7 -
                              p2**2*x0*x2*x5*x6*x7 - p1**2*x1*x2*x5*x6*x7 - 2*p1*p2*x1*x2*x5*x6*x7 - p2**2*x1*x2*x5*x6*x7 - p1**2*x0*x3*x5*x6*x7 - 2*p1*p2*x0*x3*x5*x6*x7 - p2**2*x0*x3*x5*x6*x7 -
                              p1**2*x1*x3*x5*x6*x7 - 2*p1*p2*x1*x3*x5*x6*x7 - p2**2*x1*x3*x5*x6*x7 - p1**2*x2*x3*x5*x6*x7 - p1**2*x0*x4*x5*x6*x7 - 2*p1*p2*x0*x4*x5*x6*x7 - p2**2*x0*x4*x5*x6*x7 -
                              p1**2*x1*x4*x5*x6*x7 - 2*p1*p2*x1*x4*x5*x6*x7 - p2**2*x1*x4*x5*x6*x7 - p1**2*x2*x4*x5*x6*x7 - 2*p1*p2*x2*x4*x5*x6*x7 - p2**2*x2*x4*x5*x6*x7 - p2**2*x3*x4*x5*x6*x7 -
                              p1**2*x0*x1*x2*x3*x8 - p1**2*x0*x1*x3*x4*x8 - p1**2*x0*x2*x3*x4*x8 - p1**2*x0*x1*x2*x5*x8 - 2*p1*p2*x0*x1*x2*x5*x8 - p2**2*x0*x1*x2*x5*x8 - p1**2*x0*x1*x3*x5*x8 -
                              2*p1*p2*x0*x1*x3*x5*x8 - p2**2*x0*x1*x3*x5*x8 - p1**2*x1*x2*x3*x5*x8 - p1**2*x0*x1*x4*x5*x8 - 2*p1*p2*x0*x1*x4*x5*x8 - p2**2*x0*x1*x4*x5*x8 - p1**2*x0*x2*x4*x5*x8 -
                              2*p1*p2*x0*x2*x4*x5*x8 - p2**2*x0*x2*x4*x5*x8 - p1**2*x0*x3*x4*x5*x8 - 2*p1*p2*x0*x3*x4*x5*x8 - p2**2*x0*x3*x4*x5*x8 - p1**2*x1*x3*x4*x5*x8 - p1**2*x2*x3*x4*x5*x8 -
                              p1**2*x0*x2*x3*x6*x8 - p1**2*x1*x2*x3*x6*x8 - p1**2*x0*x3*x4*x6*x8 - p1**2*x1*x3*x4*x6*x8 - p1**2*x2*x3*x4*x6*x8 - p1**2*x0*x2*x5*x6*x8 - 2*p1*p2*x0*x2*x5*x6*x8 -
                              p2**2*x0*x2*x5*x6*x8 - p1**2*x1*x2*x5*x6*x8 - 2*p1*p2*x1*x2*x5*x6*x8 - p2**2*x1*x2*x5*x6*x8 - p1**2*x0*x3*x5*x6*x8 - 2*p1*p2*x0*x3*x5*x6*x8 - p2**2*x0*x3*x5*x6*x8 -
                              p1**2*x1*x3*x5*x6*x8 - 2*p1*p2*x1*x3*x5*x6*x8 - p2**2*x1*x3*x5*x6*x8 - p1**2*x2*x3*x5*x6*x8 - p1**2*x0*x4*x5*x6*x8 - 2*p1*p2*x0*x4*x5*x6*x8 - p2**2*x0*x4*x5*x6*x8 -
                              p1**2*x1*x4*x5*x6*x8 - 2*p1*p2*x1*x4*x5*x6*x8 - p2**2*x1*x4*x5*x6*x8 - p1**2*x2*x4*x5*x6*x8 - 2*p1*p2*x2*x4*x5*x6*x8 - p2**2*x2*x4*x5*x6*x8 - p2**2*x3*x4*x5*x6*x8 -
                              p1**2*x0*x1*x2*x7*x8 - 2*p1*p2*x0*x1*x2*x7*x8 - p2**2*x0*x1*x2*x7*x8 - p1**2*x0*x1*x3*x7*x8 - 2*p1*p2*x0*x1*x3*x7*x8 - p2**2*x0*x1*x3*x7*x8 - p1**2*x1*x2*x3*x7*x8 -
                              p1**2*x0*x1*x4*x7*x8 - 2*p1*p2*x0*x1*x4*x7*x8 - p2**2*x0*x1*x4*x7*x8 - p1**2*x0*x2*x4*x7*x8 - 2*p1*p2*x0*x2*x4*x7*x8 - p2**2*x0*x2*x4*x7*x8 - p1**2*x0*x3*x4*x7*x8 -
                              2*p1*p2*x0*x3*x4*x7*x8 - p2**2*x0*x3*x4*x7*x8 - p1**2*x1*x3*x4*x7*x8 - p1**2*x2*x3*x4*x7*x8 - p1**2*x0*x2*x6*x7*x8 - 2*p1*p2*x0*x2*x6*x7*x8 - p2**2*x0*x2*x6*x7*x8 -
                              p1**2*x1*x2*x6*x7*x8 - 2*p1*p2*x1*x2*x6*x7*x8 - p2**2*x1*x2*x6*x7*x8 - p1**2*x0*x3*x6*x7*x8 - 2*p1*p2*x0*x3*x6*x7*x8 - p2**2*x0*x3*x6*x7*x8 - p1**2*x1*x3*x6*x7*x8 -
                              2*p1*p2*x1*x3*x6*x7*x8 - p2**2*x1*x3*x6*x7*x8 - p1**2*x2*x3*x6*x7*x8 - p1**2*x0*x4*x6*x7*x8 - 2*p1*p2*x0*x4*x6*x7*x8 - p2**2*x0*x4*x6*x7*x8 - p1**2*x1*x4*x6*x7*x8 -
                              2*p1*p2*x1*x4*x6*x7*x8 - p2**2*x1*x4*x6*x7*x8 - p1**2*x2*x4*x6*x7*x8 - 2*p1*p2*x2*x4*x6*x7*x8 - p2**2*x2*x4*x6*x7*x8 - p2**2*x3*x4*x6*x7*x8 - p1**2*x0*x1*x2*x3*x9 -
                              p1**2*x0*x1*x3*x4*x9 - p1**2*x0*x2*x3*x4*x9 - p2**2*x0*x1*x2*x5*x9 - p1**2*x0*x1*x3*x5*x9 - 2*p1*p2*x0*x1*x3*x5*x9 - p2**2*x0*x1*x3*x5*x9 - p1**2*x1*x2*x3*x5*x9 -
                              p2**2*x0*x1*x4*x5*x9 - p2**2*x0*x2*x4*x5*x9 - p1**2*x0*x3*x4*x5*x9 - 2*p1*p2*x0*x3*x4*x5*x9 - p2**2*x0*x3*x4*x5*x9 - p1**2*x1*x3*x4*x5*x9 - p1**2*x2*x3*x4*x5*x9 -
                              p1**2*x0*x2*x3*x6*x9 - p1**2*x1*x2*x3*x6*x9 - p1**2*x0*x3*x4*x6*x9 - p1**2*x1*x3*x4*x6*x9 - p1**2*x2*x3*x4*x6*x9 - p2**2*x0*x2*x5*x6*x9 - p2**2*x1*x2*x5*x6*x9 -
                              p1**2*x0*x3*x5*x6*x9 - 2*p1*p2*x0*x3*x5*x6*x9 - p2**2*x0*x3*x5*x6*x9 - p1**2*x1*x3*x5*x6*x9 - 2*p1*p2*x1*x3*x5*x6*x9 - p2**2*x1*x3*x5*x6*x9 - p1**2*x2*x3*x5*x6*x9 -
                              p2**2*x0*x4*x5*x6*x9 - p2**2*x1*x4*x5*x6*x9 - p2**2*x2*x4*x5*x6*x9 - p2**2*x3*x4*x5*x6*x9 - p1**2*x0*x1*x2*x7*x9 - 2*p1*p2*x0*x1*x2*x7*x9 - p2**2*x0*x1*x2*x7*x9 -
                              p1**2*x0*x1*x3*x7*x9 - 2*p1*p2*x0*x1*x3*x7*x9 - p2**2*x0*x1*x3*x7*x9 - p1**2*x1*x2*x3*x7*x9 - p1**2*x0*x1*x4*x7*x9 - 2*p1*p2*x0*x1*x4*x7*x9 - p2**2*x0*x1*x4*x7*x9 -
                              p1**2*x0*x2*x4*x7*x9 - 2*p1*p2*x0*x2*x4*x7*x9 - p2**2*x0*x2*x4*x7*x9 - p1**2*x0*x3*x4*x7*x9 - 2*p1*p2*x0*x3*x4*x7*x9 - p2**2*x0*x3*x4*x7*x9 - p1**2*x1*x3*x4*x7*x9 -
                              p1**2*x2*x3*x4*x7*x9 - p1**2*x0*x1*x5*x7*x9 - 2*p1*p2*x0*x1*x5*x7*x9 - p2**2*x0*x1*x5*x7*x9 - p1**2*x1*x2*x5*x7*x9 - p1**2*x0*x4*x5*x7*x9 - 2*p1*p2*x0*x4*x5*x7*x9 -
                              p2**2*x0*x4*x5*x7*x9 - p1**2*x1*x4*x5*x7*x9 - p1**2*x2*x4*x5*x7*x9 - p1**2*x0*x2*x6*x7*x9 - 2*p1*p2*x0*x2*x6*x7*x9 - p2**2*x0*x2*x6*x7*x9 - p1**2*x1*x2*x6*x7*x9 -
                              2*p1*p2*x1*x2*x6*x7*x9 - p2**2*x1*x2*x6*x7*x9 - p1**2*x0*x3*x6*x7*x9 - 2*p1*p2*x0*x3*x6*x7*x9 - p2**2*x0*x3*x6*x7*x9 - p1**2*x1*x3*x6*x7*x9 -
                              2*p1*p2*x1*x3*x6*x7*x9 - p2**2*x1*x3*x6*x7*x9 - p1**2*x2*x3*x6*x7*x9 - p1**2*x0*x4*x6*x7*x9 - 2*p1*p2*x0*x4*x6*x7*x9 - p2**2*x0*x4*x6*x7*x9 - p1**2*x1*x4*x6*x7*x9 -
                              2*p1*p2*x1*x4*x6*x7*x9 - p2**2*x1*x4*x6*x7*x9 - p1**2*x2*x4*x6*x7*x9 - 2*p1*p2*x2*x4*x6*x7*x9 - p2**2*x2*x4*x6*x7*x9 - p2**2*x3*x4*x6*x7*x9 - p1**2*x0*x5*x6*x7*x9 -
                              2*p1*p2*x0*x5*x6*x7*x9 - p2**2*x0*x5*x6*x7*x9 - p1**2*x1*x5*x6*x7*x9 - 2*p1*p2*x1*x5*x6*x7*x9 - p2**2*x1*x5*x6*x7*x9 - p1**2*x2*x5*x6*x7*x9 - p2**2*x4*x5*x6*x7*x9 -
                              p1**2*x0*x1*x2*x8*x9 - p1**2*x0*x1*x4*x8*x9 - p1**2*x0*x2*x4*x8*x9 - p1**2*x0*x1*x5*x8*x9 - 2*p1*p2*x0*x1*x5*x8*x9 - p2**2*x0*x1*x5*x8*x9 - p1**2*x1*x2*x5*x8*x9 -
                              p1**2*x0*x4*x5*x8*x9 - 2*p1*p2*x0*x4*x5*x8*x9 - p2**2*x0*x4*x5*x8*x9 - p1**2*x1*x4*x5*x8*x9 - p1**2*x2*x4*x5*x8*x9 - p1**2*x0*x2*x6*x8*x9 - p1**2*x1*x2*x6*x8*x9 -
                              p1**2*x0*x4*x6*x8*x9 - p1**2*x1*x4*x6*x8*x9 - p1**2*x2*x4*x6*x8*x9 - p1**2*x0*x5*x6*x8*x9 - 2*p1*p2*x0*x5*x6*x8*x9 - p2**2*x0*x5*x6*x8*x9 - p1**2*x1*x5*x6*x8*x9 -
                              2*p1*p2*x1*x5*x6*x8*x9 - p2**2*x1*x5*x6*x8*x9 - p1**2*x2*x5*x6*x8*x9 - p2**2*x4*x5*x6*x8*x9 - p1**2*x0*x1*x7*x8*x9 - 2*p1*p2*x0*x1*x7*x8*x9 - p2**2*x0*x1*x7*x8*x9 -
                              p1**2*x1*x2*x7*x8*x9 - p1**2*x0*x4*x7*x8*x9 - 2*p1*p2*x0*x4*x7*x8*x9 - p2**2*x0*x4*x7*x8*x9 - p1**2*x1*x4*x7*x8*x9 - p1**2*x2*x4*x7*x8*x9 - p1**2*x0*x6*x7*x8*x9 -
                              2*p1*p2*x0*x6*x7*x8*x9 - p2**2*x0*x6*x7*x8*x9 - p1**2*x1*x6*x7*x8*x9 - 2*p1*p2*x1*x6*x7*x8*x9 - p2**2*x1*x6*x7*x8*x9 - p1**2*x2*x6*x7*x8*x9 - p2**2*x4*x6*x7*x8*x9 -
                              p1**2*x0*x1*x2*x3*x10 - p1**2*x0*x1*x3*x4*x10 - p1**2*x0*x2*x3*x4*x10 - p2**2*x0*x1*x2*x5*x10 - p1**2*x0*x1*x3*x5*x10 - 2*p1*p2*x0*x1*x3*x5*x10 -
                              p2**2*x0*x1*x3*x5*x10 - p1**2*x0*x2*x3*x5*x10 - p2**2*x0*x1*x4*x5*x10 - p2**2*x0*x2*x4*x5*x10 - p2**2*x0*x3*x4*x5*x10 - p1**2*x0*x2*x3*x6*x10 - p1**2*x1*x2*x3*x6*x10 -
                              p1**2*x0*x3*x4*x6*x10 - p1**2*x1*x3*x4*x6*x10 - p1**2*x2*x3*x4*x6*x10 - p2**2*x0*x2*x5*x6*x10 - p2**2*x1*x2*x5*x6*x10 - p1**2*x0*x3*x5*x6*x10 -
                              2*p1*p2*x0*x3*x5*x6*x10 - p2**2*x0*x3*x5*x6*x10 - p1**2*x1*x3*x5*x6*x10 - 2*p1*p2*x1*x3*x5*x6*x10 - p2**2*x1*x3*x5*x6*x10 - p1**2*x2*x3*x5*x6*x10 -
                              p2**2*x0*x4*x5*x6*x10 - p2**2*x1*x4*x5*x6*x10 - p2**2*x2*x4*x5*x6*x10 - p2**2*x3*x4*x5*x6*x10 - p1**2*x0*x1*x2*x7*x10 - 2*p1*p2*x0*x1*x2*x7*x10 -
                              p2**2*x0*x1*x2*x7*x10 - p1**2*x0*x1*x3*x7*x10 - 2*p1*p2*x0*x1*x3*x7*x10 - p2**2*x0*x1*x3*x7*x10 - p1**2*x1*x2*x3*x7*x10 - p1**2*x0*x1*x4*x7*x10 -
                              2*p1*p2*x0*x1*x4*x7*x10 - p2**2*x0*x1*x4*x7*x10 - p1**2*x0*x2*x4*x7*x10 - 2*p1*p2*x0*x2*x4*x7*x10 - p2**2*x0*x2*x4*x7*x10 - p1**2*x0*x3*x4*x7*x10 -
                              2*p1*p2*x0*x3*x4*x7*x10 - p2**2*x0*x3*x4*x7*x10 - p1**2*x1*x3*x4*x7*x10 - p1**2*x2*x3*x4*x7*x10 - p1**2*x0*x1*x5*x7*x10 - 2*p1*p2*x0*x1*x5*x7*x10 -
                              p2**2*x0*x1*x5*x7*x10 - p1**2*x0*x2*x5*x7*x10 - 2*p1*p2*x0*x2*x5*x7*x10 - p2**2*x0*x2*x5*x7*x10 - p2**2*x1*x2*x5*x7*x10 - p1**2*x0*x3*x5*x7*x10 -
                              2*p1*p2*x0*x3*x5*x7*x10 - p2**2*x0*x3*x5*x7*x10 - p1**2*x1*x3*x5*x7*x10 - 2*p1*p2*x1*x3*x5*x7*x10 - p2**2*x1*x3*x5*x7*x10 - p1**2*x2*x3*x5*x7*x10 -
                              p2**2*x1*x4*x5*x7*x10 - p2**2*x2*x4*x5*x7*x10 - p2**2*x3*x4*x5*x7*x10 - p1**2*x0*x2*x6*x7*x10 - 2*p1*p2*x0*x2*x6*x7*x10 - p2**2*x0*x2*x6*x7*x10 -
                              p1**2*x1*x2*x6*x7*x10 - 2*p1*p2*x1*x2*x6*x7*x10 - p2**2*x1*x2*x6*x7*x10 - p1**2*x0*x3*x6*x7*x10 - 2*p1*p2*x0*x3*x6*x7*x10 - p2**2*x0*x3*x6*x7*x10 -
                              p1**2*x1*x3*x6*x7*x10 - 2*p1*p2*x1*x3*x6*x7*x10 - p2**2*x1*x3*x6*x7*x10 - p1**2*x2*x3*x6*x7*x10 - p1**2*x0*x4*x6*x7*x10 - 2*p1*p2*x0*x4*x6*x7*x10 -
                              p2**2*x0*x4*x6*x7*x10 - p1**2*x1*x4*x6*x7*x10 - 2*p1*p2*x1*x4*x6*x7*x10 - p2**2*x1*x4*x6*x7*x10 - p1**2*x2*x4*x6*x7*x10 - 2*p1*p2*x2*x4*x6*x7*x10 -
                              p2**2*x2*x4*x6*x7*x10 - p2**2*x3*x4*x6*x7*x10 - p1**2*x0*x5*x6*x7*x10 - 2*p1*p2*x0*x5*x6*x7*x10 - p2**2*x0*x5*x6*x7*x10 - p1**2*x1*x5*x6*x7*x10 -
                              2*p1*p2*x1*x5*x6*x7*x10 - p2**2*x1*x5*x6*x7*x10 - p1**2*x2*x5*x6*x7*x10 - p2**2*x4*x5*x6*x7*x10 - p1**2*x0*x1*x2*x8*x10 - p1**2*x0*x2*x3*x8*x10 -
                              p1**2*x1*x2*x3*x8*x10 - p1**2*x0*x1*x4*x8*x10 - p1**2*x0*x2*x4*x8*x10 - p1**2*x0*x3*x4*x8*x10 - p1**2*x1*x3*x4*x8*x10 - p1**2*x2*x3*x4*x8*x10 - p1**2*x0*x1*x5*x8*x10 -
                              2*p1*p2*x0*x1*x5*x8*x10 - p2**2*x0*x1*x5*x8*x10 - p1**2*x0*x2*x5*x8*x10 - 2*p1*p2*x0*x2*x5*x8*x10 - p2**2*x0*x2*x5*x8*x10 - p2**2*x1*x2*x5*x8*x10 -
                              p1**2*x0*x3*x5*x8*x10 - 2*p1*p2*x0*x3*x5*x8*x10 - p2**2*x0*x3*x5*x8*x10 - p1**2*x1*x3*x5*x8*x10 - 2*p1*p2*x1*x3*x5*x8*x10 - p2**2*x1*x3*x5*x8*x10 -
                              p1**2*x2*x3*x5*x8*x10 - p2**2*x1*x4*x5*x8*x10 - p2**2*x2*x4*x5*x8*x10 - p2**2*x3*x4*x5*x8*x10 - p1**2*x0*x2*x6*x8*x10 - p1**2*x1*x2*x6*x8*x10 - p1**2*x0*x4*x6*x8*x10 -
                              p1**2*x1*x4*x6*x8*x10 - p1**2*x2*x4*x6*x8*x10 - p1**2*x0*x5*x6*x8*x10 - 2*p1*p2*x0*x5*x6*x8*x10 - p2**2*x0*x5*x6*x8*x10 - p1**2*x1*x5*x6*x8*x10 -
                              2*p1*p2*x1*x5*x6*x8*x10 - p2**2*x1*x5*x6*x8*x10 - p1**2*x2*x5*x6*x8*x10 - p2**2*x4*x5*x6*x8*x10 - p1**2*x0*x1*x7*x8*x10 - 2*p1*p2*x0*x1*x7*x8*x10 -
                              p2**2*x0*x1*x7*x8*x10 - p1**2*x0*x2*x7*x8*x10 - 2*p1*p2*x0*x2*x7*x8*x10 - p2**2*x0*x2*x7*x8*x10 - p2**2*x1*x2*x7*x8*x10 - p1**2*x0*x3*x7*x8*x10 -
                              2*p1*p2*x0*x3*x7*x8*x10 - p2**2*x0*x3*x7*x8*x10 - p1**2*x1*x3*x7*x8*x10 - 2*p1*p2*x1*x3*x7*x8*x10 - p2**2*x1*x3*x7*x8*x10 - p1**2*x2*x3*x7*x8*x10 -
                              p2**2*x1*x4*x7*x8*x10 - p2**2*x2*x4*x7*x8*x10 - p2**2*x3*x4*x7*x8*x10 - p1**2*x0*x6*x7*x8*x10 - 2*p1*p2*x0*x6*x7*x8*x10 - p2**2*x0*x6*x7*x8*x10 -
                              p1**2*x1*x6*x7*x8*x10 - 2*p1*p2*x1*x6*x7*x8*x10 - p2**2*x1*x6*x7*x8*x10 - p1**2*x2*x6*x7*x8*x10 - p2**2*x4*x6*x7*x8*x10 - p1**2*x0*x2*x3*x9*x10 -
                              p1**2*x1*x2*x3*x9*x10 - p1**2*x0*x3*x4*x9*x10 - p1**2*x1*x3*x4*x9*x10 - p1**2*x2*x3*x4*x9*x10 - p2**2*x0*x2*x5*x9*x10 - p2**2*x1*x2*x5*x9*x10 - p1**2*x0*x3*x5*x9*x10 -
                              2*p1*p2*x0*x3*x5*x9*x10 - p2**2*x0*x3*x5*x9*x10 - p1**2*x1*x3*x5*x9*x10 - 2*p1*p2*x1*x3*x5*x9*x10 - p2**2*x1*x3*x5*x9*x10 - p1**2*x2*x3*x5*x9*x10 -
                              p2**2*x0*x4*x5*x9*x10 - p2**2*x1*x4*x5*x9*x10 - p2**2*x2*x4*x5*x9*x10 - p2**2*x3*x4*x5*x9*x10 - p1**2*x0*x2*x7*x9*x10 - 2*p1*p2*x0*x2*x7*x9*x10 -
                              p2**2*x0*x2*x7*x9*x10 - p1**2*x1*x2*x7*x9*x10 - 2*p1*p2*x1*x2*x7*x9*x10 - p2**2*x1*x2*x7*x9*x10 - p1**2*x0*x3*x7*x9*x10 - 2*p1*p2*x0*x3*x7*x9*x10 -
                              p2**2*x0*x3*x7*x9*x10 - p1**2*x1*x3*x7*x9*x10 - 2*p1*p2*x1*x3*x7*x9*x10 - p2**2*x1*x3*x7*x9*x10 - p1**2*x2*x3*x7*x9*x10 - p1**2*x0*x4*x7*x9*x10 -
                              2*p1*p2*x0*x4*x7*x9*x10 - p2**2*x0*x4*x7*x9*x10 - p1**2*x1*x4*x7*x9*x10 - 2*p1*p2*x1*x4*x7*x9*x10 - p2**2*x1*x4*x7*x9*x10 - p1**2*x2*x4*x7*x9*x10 -
                              2*p1*p2*x2*x4*x7*x9*x10 - p2**2*x2*x4*x7*x9*x10 - p2**2*x3*x4*x7*x9*x10 - p1**2*x0*x5*x7*x9*x10 - 2*p1*p2*x0*x5*x7*x9*x10 - p2**2*x0*x5*x7*x9*x10 -
                              p1**2*x1*x5*x7*x9*x10 - 2*p1*p2*x1*x5*x7*x9*x10 - p2**2*x1*x5*x7*x9*x10 - p1**2*x2*x5*x7*x9*x10 - p2**2*x4*x5*x7*x9*x10 - p1**2*x0*x2*x8*x9*x10 -
                              p1**2*x1*x2*x8*x9*x10 - p1**2*x0*x4*x8*x9*x10 - p1**2*x1*x4*x8*x9*x10 - p1**2*x2*x4*x8*x9*x10 - p1**2*x0*x5*x8*x9*x10 - 2*p1*p2*x0*x5*x8*x9*x10 -
                              p2**2*x0*x5*x8*x9*x10 - p1**2*x1*x5*x8*x9*x10 - 2*p1*p2*x1*x5*x8*x9*x10 - p2**2*x1*x5*x8*x9*x10 - p1**2*x2*x5*x8*x9*x10 - p2**2*x4*x5*x8*x9*x10 -
                              p1**2*x0*x7*x8*x9*x10 - 2*p1*p2*x0*x7*x8*x9*x10 - p2**2*x0*x7*x8*x9*x10 - p1**2*x1*x7*x8*x9*x10 - 2*p1*p2*x1*x7*x8*x9*x10 - p2**2*x1*x7*x8*x9*x10 -
                              p1**2*x2*x7*x8*x9*x10 - p2**2*x4*x7*x8*x9*x10 - p1**2*x0*x1*x2*x3*x11 - p1**2*x0*x1*x3*x4*x11 - p1**2*x0*x2*x3*x4*x11 - p1**2*x0*x1*x3*x5*x11 - p1**2*x0*x2*x3*x5*x11 -
                              p1**2*x0*x2*x3*x6*x11 - p1**2*x1*x2*x3*x6*x11 - p1**2*x0*x3*x4*x6*x11 - p1**2*x1*x3*x4*x6*x11 - p1**2*x2*x3*x4*x6*x11 - p1**2*x0*x3*x5*x6*x11 - p1**2*x1*x3*x5*x6*x11 -
                              p1**2*x2*x3*x5*x6*x11 - p1**2*x0*x1*x2*x7*x11 - 2*p1*p2*x0*x1*x2*x7*x11 - p2**2*x0*x1*x2*x7*x11 - p1**2*x0*x1*x3*x7*x11 - 2*p1*p2*x0*x1*x3*x7*x11 -
                              p2**2*x0*x1*x3*x7*x11 - p1**2*x1*x2*x3*x7*x11 - p1**2*x0*x1*x4*x7*x11 - 2*p1*p2*x0*x1*x4*x7*x11 - p2**2*x0*x1*x4*x7*x11 - p1**2*x0*x2*x4*x7*x11 -
                              2*p1*p2*x0*x2*x4*x7*x11 - p2**2*x0*x2*x4*x7*x11 - p1**2*x0*x3*x4*x7*x11 - 2*p1*p2*x0*x3*x4*x7*x11 - p2**2*x0*x3*x4*x7*x11 - p1**2*x1*x3*x4*x7*x11 -
                              p1**2*x2*x3*x4*x7*x11 - p1**2*x0*x1*x5*x7*x11 - 2*p1*p2*x0*x1*x5*x7*x11 - p2**2*x0*x1*x5*x7*x11 - p1**2*x0*x2*x5*x7*x11 - 2*p1*p2*x0*x2*x5*x7*x11 -
                              p2**2*x0*x2*x5*x7*x11 - p1**2*x0*x3*x5*x7*x11 - 2*p1*p2*x0*x3*x5*x7*x11 - p2**2*x0*x3*x5*x7*x11 - p1**2*x1*x3*x5*x7*x11 - p1**2*x2*x3*x5*x7*x11 -
                              p1**2*x0*x2*x6*x7*x11 - 2*p1*p2*x0*x2*x6*x7*x11 - p2**2*x0*x2*x6*x7*x11 - p1**2*x1*x2*x6*x7*x11 - 2*p1*p2*x1*x2*x6*x7*x11 - p2**2*x1*x2*x6*x7*x11 -
                              p1**2*x0*x3*x6*x7*x11 - 2*p1*p2*x0*x3*x6*x7*x11 - p2**2*x0*x3*x6*x7*x11 - p1**2*x1*x3*x6*x7*x11 - 2*p1*p2*x1*x3*x6*x7*x11 - p2**2*x1*x3*x6*x7*x11 -
                              p1**2*x2*x3*x6*x7*x11 - p1**2*x0*x4*x6*x7*x11 - 2*p1*p2*x0*x4*x6*x7*x11 - p2**2*x0*x4*x6*x7*x11 - p1**2*x1*x4*x6*x7*x11 - 2*p1*p2*x1*x4*x6*x7*x11 -
                              p2**2*x1*x4*x6*x7*x11 - p1**2*x2*x4*x6*x7*x11 - 2*p1*p2*x2*x4*x6*x7*x11 - p2**2*x2*x4*x6*x7*x11 - p2**2*x3*x4*x6*x7*x11 - p1**2*x0*x5*x6*x7*x11 -
                              2*p1*p2*x0*x5*x6*x7*x11 - p2**2*x0*x5*x6*x7*x11 - p1**2*x1*x5*x6*x7*x11 - 2*p1*p2*x1*x5*x6*x7*x11 - p2**2*x1*x5*x6*x7*x11 - p1**2*x2*x5*x6*x7*x11 -
                              2*p1*p2*x2*x5*x6*x7*x11 - p2**2*x2*x5*x6*x7*x11 - p2**2*x3*x5*x6*x7*x11 - p1**2*x0*x1*x2*x8*x11 - 2*p1*p2*x0*x1*x2*x8*x11 - p2**2*x0*x1*x2*x8*x11 -
                              p2**2*x0*x1*x3*x8*x11 - p1**2*x0*x2*x3*x8*x11 - p1**2*x1*x2*x3*x8*x11 - p1**2*x0*x1*x4*x8*x11 - 2*p1*p2*x0*x1*x4*x8*x11 - p2**2*x0*x1*x4*x8*x11 -
                              p1**2*x0*x2*x4*x8*x11 - 2*p1*p2*x0*x2*x4*x8*x11 - p2**2*x0*x2*x4*x8*x11 - p1**2*x0*x3*x4*x8*x11 - 2*p1*p2*x0*x3*x4*x8*x11 - p2**2*x0*x3*x4*x8*x11 -
                              p1**2*x1*x3*x4*x8*x11 - p1**2*x2*x3*x4*x8*x11 - p1**2*x0*x1*x5*x8*x11 - 2*p1*p2*x0*x1*x5*x8*x11 - p2**2*x0*x1*x5*x8*x11 - p1**2*x0*x2*x5*x8*x11 -
                              2*p1*p2*x0*x2*x5*x8*x11 - p2**2*x0*x2*x5*x8*x11 - p1**2*x0*x3*x5*x8*x11 - 2*p1*p2*x0*x3*x5*x8*x11 - p2**2*x0*x3*x5*x8*x11 - p1**2*x1*x3*x5*x8*x11 -
                              p1**2*x2*x3*x5*x8*x11 - p1**2*x0*x2*x6*x8*x11 - 2*p1*p2*x0*x2*x6*x8*x11 - p2**2*x0*x2*x6*x8*x11 - p1**2*x1*x2*x6*x8*x11 - 2*p1*p2*x1*x2*x6*x8*x11 -
                              p2**2*x1*x2*x6*x8*x11 - p2**2*x0*x3*x6*x8*x11 - p2**2*x1*x3*x6*x8*x11 - p1**2*x0*x4*x6*x8*x11 - 2*p1*p2*x0*x4*x6*x8*x11 - p2**2*x0*x4*x6*x8*x11 -
                              p1**2*x1*x4*x6*x8*x11 - 2*p1*p2*x1*x4*x6*x8*x11 - p2**2*x1*x4*x6*x8*x11 - p1**2*x2*x4*x6*x8*x11 - 2*p1*p2*x2*x4*x6*x8*x11 - p2**2*x2*x4*x6*x8*x11 -
                              p2**2*x3*x4*x6*x8*x11 - p1**2*x0*x5*x6*x8*x11 - 2*p1*p2*x0*x5*x6*x8*x11 - p2**2*x0*x5*x6*x8*x11 - p1**2*x1*x5*x6*x8*x11 - 2*p1*p2*x1*x5*x6*x8*x11 -
                              p2**2*x1*x5*x6*x8*x11 - p1**2*x2*x5*x6*x8*x11 - 2*p1*p2*x2*x5*x6*x8*x11 - p2**2*x2*x5*x6*x8*x11 - p2**2*x3*x5*x6*x8*x11 - p1**2*x0*x1*x7*x8*x11 -
                              2*p1*p2*x0*x1*x7*x8*x11 - p2**2*x0*x1*x7*x8*x11 - p1**2*x0*x2*x7*x8*x11 - 2*p1*p2*x0*x2*x7*x8*x11 - p2**2*x0*x2*x7*x8*x11 - p1**2*x0*x3*x7*x8*x11 -
                              2*p1*p2*x0*x3*x7*x8*x11 - p2**2*x0*x3*x7*x8*x11 - p1**2*x1*x3*x7*x8*x11 - p1**2*x2*x3*x7*x8*x11 - p1**2*x0*x6*x7*x8*x11 - 2*p1*p2*x0*x6*x7*x8*x11 -
                              p2**2*x0*x6*x7*x8*x11 - p1**2*x1*x6*x7*x8*x11 - 2*p1*p2*x1*x6*x7*x8*x11 - p2**2*x1*x6*x7*x8*x11 - p1**2*x2*x6*x7*x8*x11 - 2*p1*p2*x2*x6*x7*x8*x11 -
                              p2**2*x2*x6*x7*x8*x11 - p2**2*x3*x6*x7*x8*x11 - p2**2*x0*x1*x2*x9*x11 - p2**2*x0*x1*x3*x9*x11 - p1**2*x0*x2*x3*x9*x11 - p1**2*x1*x2*x3*x9*x11 - p2**2*x0*x1*x4*x9*x11 -
                              p2**2*x0*x2*x4*x9*x11 - p1**2*x0*x3*x4*x9*x11 - 2*p1*p2*x0*x3*x4*x9*x11 - p2**2*x0*x3*x4*x9*x11 - p1**2*x1*x3*x4*x9*x11 - p1**2*x2*x3*x4*x9*x11 -
                              p2**2*x0*x1*x5*x9*x11 - p2**2*x0*x2*x5*x9*x11 - p1**2*x0*x3*x5*x9*x11 - 2*p1*p2*x0*x3*x5*x9*x11 - p2**2*x0*x3*x5*x9*x11 - p1**2*x1*x3*x5*x9*x11 -
                              p1**2*x2*x3*x5*x9*x11 - p2**2*x0*x2*x6*x9*x11 - p2**2*x1*x2*x6*x9*x11 - p2**2*x0*x3*x6*x9*x11 - p2**2*x1*x3*x6*x9*x11 - p2**2*x0*x4*x6*x9*x11 - p2**2*x1*x4*x6*x9*x11 -
                              p2**2*x2*x4*x6*x9*x11 - p2**2*x3*x4*x6*x9*x11 - p2**2*x0*x5*x6*x9*x11 - p2**2*x1*x5*x6*x9*x11 - p2**2*x2*x5*x6*x9*x11 - p2**2*x3*x5*x6*x9*x11 - p1**2*x0*x2*x7*x9*x11 -
                              2*p1*p2*x0*x2*x7*x9*x11 - p2**2*x0*x2*x7*x9*x11 - p1**2*x1*x2*x7*x9*x11 - p1**2*x0*x3*x7*x9*x11 - 2*p1*p2*x0*x3*x7*x9*x11 - p2**2*x0*x3*x7*x9*x11 -
                              p1**2*x1*x3*x7*x9*x11 - p1**2*x2*x3*x7*x9*x11 - p1**2*x0*x4*x7*x9*x11 - 2*p1*p2*x0*x4*x7*x9*x11 - p2**2*x0*x4*x7*x9*x11 - p1**2*x1*x4*x7*x9*x11 -
                              p1**2*x2*x4*x7*x9*x11 - p1**2*x0*x5*x7*x9*x11 - 2*p1*p2*x0*x5*x7*x9*x11 - p2**2*x0*x5*x7*x9*x11 - p1**2*x1*x5*x7*x9*x11 - p1**2*x2*x5*x7*x9*x11 -
                              p2**2*x2*x6*x7*x9*x11 - p2**2*x3*x6*x7*x9*x11 - p2**2*x4*x6*x7*x9*x11 - p2**2*x5*x6*x7*x9*x11 - p2**2*x0*x1*x8*x9*x11 - p1**2*x0*x2*x8*x9*x11 - p1**2*x1*x2*x8*x9*x11 -
                              p1**2*x0*x4*x8*x9*x11 - 2*p1*p2*x0*x4*x8*x9*x11 - p2**2*x0*x4*x8*x9*x11 - p1**2*x1*x4*x8*x9*x11 - p1**2*x2*x4*x8*x9*x11 - p1**2*x0*x5*x8*x9*x11 -
                              2*p1*p2*x0*x5*x8*x9*x11 - p2**2*x0*x5*x8*x9*x11 - p1**2*x1*x5*x8*x9*x11 - p1**2*x2*x5*x8*x9*x11 - p2**2*x0*x6*x8*x9*x11 - p2**2*x1*x6*x8*x9*x11 -
                              p2**2*x4*x6*x8*x9*x11 - p2**2*x5*x6*x8*x9*x11 - p1**2*x0*x7*x8*x9*x11 - 2*p1*p2*x0*x7*x8*x9*x11 - p2**2*x0*x7*x8*x9*x11 - p1**2*x1*x7*x8*x9*x11 -
                              p1**2*x2*x7*x8*x9*x11 - p2**2*x6*x7*x8*x9*x11 - p2**2*x0*x1*x2*x10*x11 - p2**2*x0*x1*x3*x10*x11 - p2**2*x0*x1*x4*x10*x11 - p2**2*x0*x2*x4*x10*x11 -
                              p2**2*x0*x3*x4*x10*x11 - p2**2*x0*x1*x5*x10*x11 - p2**2*x0*x2*x5*x10*x11 - p2**2*x0*x3*x5*x10*x11 - p2**2*x0*x2*x6*x10*x11 - p2**2*x1*x2*x6*x10*x11 -
                              p2**2*x0*x3*x6*x10*x11 - p2**2*x1*x3*x6*x10*x11 - p2**2*x0*x4*x6*x10*x11 - p2**2*x1*x4*x6*x10*x11 - p2**2*x2*x4*x6*x10*x11 - p2**2*x3*x4*x6*x10*x11 -
                              p2**2*x0*x5*x6*x10*x11 - p2**2*x1*x5*x6*x10*x11 - p2**2*x2*x5*x6*x10*x11 - p2**2*x3*x5*x6*x10*x11 - p2**2*x1*x2*x7*x10*x11 - p2**2*x1*x3*x7*x10*x11 -
                              p2**2*x1*x4*x7*x10*x11 - p2**2*x2*x4*x7*x10*x11 - p2**2*x3*x4*x7*x10*x11 - p2**2*x1*x5*x7*x10*x11 - p2**2*x2*x5*x7*x10*x11 - p2**2*x3*x5*x7*x10*x11 -
                              p2**2*x2*x6*x7*x10*x11 - p2**2*x3*x6*x7*x10*x11 - p2**2*x4*x6*x7*x10*x11 - p2**2*x5*x6*x7*x10*x11 - p2**2*x0*x1*x8*x10*x11 - p2**2*x0*x2*x8*x10*x11 -
                              p2**2*x1*x2*x8*x10*x11 - p2**2*x0*x3*x8*x10*x11 - p2**2*x1*x3*x8*x10*x11 - p2**2*x1*x4*x8*x10*x11 - p2**2*x2*x4*x8*x10*x11 - p2**2*x3*x4*x8*x10*x11 -
                              p2**2*x1*x5*x8*x10*x11 - p2**2*x2*x5*x8*x10*x11 - p2**2*x3*x5*x8*x10*x11 - p2**2*x0*x6*x8*x10*x11 - p2**2*x1*x6*x8*x10*x11 - p2**2*x4*x6*x8*x10*x11 -
                              p2**2*x5*x6*x8*x10*x11 - p2**2*x1*x7*x8*x10*x11 - p2**2*x2*x7*x8*x10*x11 - p2**2*x3*x7*x8*x10*x11 - p2**2*x6*x7*x8*x10*x11 - p2**2*x0*x2*x9*x10*x11 -
                              p2**2*x1*x2*x9*x10*x11 - p2**2*x0*x3*x9*x10*x11 - p2**2*x1*x3*x9*x10*x11 - p2**2*x0*x4*x9*x10*x11 - p2**2*x1*x4*x9*x10*x11 - p2**2*x2*x4*x9*x10*x11 -
                              p2**2*x3*x4*x9*x10*x11 - p2**2*x0*x5*x9*x10*x11 - p2**2*x1*x5*x9*x10*x11 - p2**2*x2*x5*x9*x10*x11 - p2**2*x3*x5*x9*x10*x11 - p2**2*x2*x7*x9*x10*x11 -
                              p2**2*x3*x7*x9*x10*x11 - p2**2*x4*x7*x9*x10*x11 - p2**2*x5*x7*x9*x10*x11 - p2**2*x0*x8*x9*x10*x11 - p2**2*x1*x8*x9*x10*x11 - p2**2*x4*x8*x9*x10*x11 -
                              p2**2*x5*x8*x9*x10*x11 - p2**2*x7*x8*x9*x10*x11""")

    #@attr('active')
    def test_error_messages(self):
        self.assertRaisesRegexp(AssertionError, '(M|m)ismatch.*propagators.*Feynman_parameters', \
                                LoopIntegralFromPropagators, ['k1**2'], ['k1'], Feynman_parameters=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*loop_momenta.*symbol', \
                                LoopIntegralFromPropagators, ['k1**2', '(k1+k2)**2'], ['k1', '2*k2'], Feynman_parameters=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*external_momenta.*symbol', \
                                LoopIntegralFromPropagators, ['k1**2'], ['k1', 'k2'], ['alpha*p1'], Feynman_parameters=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*Lorentz_indices.*symbol.*number', \
                                LoopIntegralFromPropagators, ['k1**2'], ['k1', 'k2'], Lorentz_indices=['alpha*delta'])
        self.assertRaisesRegexp(AssertionError, '.*metric_tensor.*symbol', \
                                LoopIntegralFromPropagators, ['k1**2'], ['k1', 'k2'], metric_tensor='4')
        self.assertRaisesRegexp(AssertionError, '.*propagators.*at most quadratic', \
                                LoopIntegralFromPropagators, ['k1**2', 'p1*(k1+k2)**2'], ['k1','k2'], ['p1','p2'])
        self.assertRaisesRegexp(AssertionError, '.*replacement_rules.*list of pairs', \
                                LoopIntegralFromPropagators, ['k1**2', '(k1+k2)**2'], ['k1','k2'], replacement_rules=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*replacement_rules.*list of pairs', \
                                LoopIntegralFromPropagators, ['k1**2', '(k1+k2)**2'], ['k1','k2'], replacement_rules=[('k1*k1', 0,'z1')])
        self.assertRaisesRegexp(AssertionError, '.*replacement_rules.*at most quadratic', \
                                LoopIntegralFromPropagators, ['k1**2', '(k1+k2)**2'], ['k1','k2'], replacement_rules=[('k1*k1*k2', 0)])

    #@attr('active')
    def test_replacement_rules_box_2L(self):
        k1, k2, p1, p2, p3, p4, s, t = sp.symbols('k1 k2 p1 p2 p3 p4 s t')

        loop_momenta = [k1, k2]
        propagators = [k1**2,(k1+p2)**2,(k1-p1)**2,(k1-k2)**2,(k2+p2)**2,(k2-p1)**2,(k2+p2+p3)**2]
        replacement_rules = [( p1*p1 , 0),
                             ('p2*p2', 0),
                             ( p3*p3 , 0),
                             ('p4*p4', 0),
                             ( p1*p2 , s/2),
                             ('p2*p3', 't/2'),
                             ( p1*p3 , -s/2-t/2)]
        loop_integral = LoopIntegralFromPropagators(propagators, loop_momenta, replacement_rules=replacement_rules, Feynman_parameters='z')
        U = sympify_expression(loop_integral.U)
        F = sympify_expression(loop_integral.F)

        target_U = sympify_expression('z3*(z4 + z5 + z6) + z0*(z3 + z4 + z5 + z6) + z1*(z3 + z4 + z5 + z6) + z2*(z3 + z4 + z5 + z6)')
        target_F = sympify_expression('-(s*z3*z4*z5) + s*z2*(-(z3*z4) - z4*z5) - z0*(s*z4*z5 + t*z3*z6) - s*z1*((z3 + z4)*z5 + z2*(z3 + z4 + z5 + z6))')

        self.assertEqual( (target_U - U).simplify() , 0)
        self.assertEqual( (target_F - F).simplify() , 0)

        L = 2
        D = sympify_expression('4-2*eps')
        N_nu = 7
        target_exponent_U = N_nu - D/2 * (L+1)
        target_exponent_F = -(N_nu - D/2 * L)
        self.assertEqual(target_exponent_U, loop_integral.exponent_U)
        self.assertEqual(target_exponent_F, loop_integral.exponent_F)

    #@attr('active')
    def test_replacement_rules_box_1L(self):
        loop_momenta = ['k1']
        props = ['(k1+p1)**2 - m**2', '(k1+p1+p2)**2', '(k1+p1+p2+p3)**2', '(k1+p1+p2+p3+p4)**2']
        replacement_rules = [('p1*p1',      's1'     ),
                             ('p2*p2',       0       ),
                             ('p3*p3',       0       ),
                             ('p3*p2',     't/2'     ),
                             ('p1*p3',   '-t/2-s/2'  ),
                             ('p1*p2',   's/2-s1/2'  ),
                             ('p4*p4',       0       ),
                             ('p1*p4',   't/2-s1/2'  ),
                             ('p2*p4', 's1/2-t/2-s/2'),
                             ('p3*p4',     's/2'     )]
        loop_integral = LoopIntegralFromPropagators(props, loop_momenta, replacement_rules=replacement_rules, Feynman_parameters=['z1','z2','z3','z4'])
        U = sympify_expression(loop_integral.U)
        F = sympify_expression(loop_integral.F)

        target_U = sympify_expression('z1 + z2 + z3 + z4')
        target_F = sympify_expression('-(t*z1*z3) - (s1*z1 + s*z2)*z4 + m**2*z1*(z1 + z2 + z3 + z4)')

        self.assertEqual( (target_U - U).simplify() , 0)
        self.assertEqual( (target_F - F).simplify() , 0)

    #@attr('active')
    def test_linear_propagator(self):
        # SecDec3 -> loop/demos/8_linearprop_1L
        loop_momenta = ['k']
        props = ['k**2','(-k+p)**2','2*k*v']
        rules = [('p*p', 'ssp1'),
                 ('v*v', 'ssp2'),
                 ('p*v', 0)]
        li = LoopIntegralFromPropagators(props, loop_momenta, replacement_rules=rules, \
                                         Feynman_parameters=['z1','z2','z3'])
        U = sympify_expression(li.U)
        F = sympify_expression(li.F)

        target_U = sympify_expression('z1 + z2')
        target_F = sympify_expression('-(ssp1*z1*z2) + ssp2*z3**2')


class TestNumerator(unittest.TestCase):
    #@attr('active')
    def test_double_index_notation(self):
        mu, nu = sp.symbols('mu nu')
        numerator = 'k1(mu)*(k2(mu) + p1(mu)) + k2(mu)*k2(nu)*p1(mu)*p2(nu)'
        indices = ['mu', 'nu']
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator, Lorentz_indices=indices)
        loop_tensors = li.numerator_loop_tensors
        external_tensors = li.numerator_external_tensors

        # Note: `sympy` may reorder the terms
        target_loop_tensors = [
                                 [(0,mu),(1,mu)], # (loop_momenta[0], index), (loop_momenta[1], index)
                                 [(0,mu)],        # (loop_momenta[0], index)
                                 [(1,mu),(1,nu)], # (loop_momenta[0], index), (loop_momenta[1], index)
                              ]
        target_external_tensors =  [
                                      [],
                                      [(0,mu)],
                                      [(0,mu),(1,nu)],
                                   ]

        self.assertEqual(len(loop_tensors), 3)
        for i in range(3):
            print(i)
            self.assertEqual(loop_tensors[i], target_loop_tensors[i])

        self.assertEqual(len(external_tensors), 3)
        for i in range(3):
            print(i)
            self.assertEqual(external_tensors[i], target_external_tensors[i])

    #@attr('active')
    def test_double_index_notation_special_case(self):
        mu = sp.symbols('mu')
        numerator = 'k1(mu)'
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator, Lorentz_indices=[mu])
        tensors = li.numerator_loop_tensors

        target_tensors = [  [(0,mu)]  ]
        self.assertEqual(tensors, target_tensors)

    def test_double_index_notation_2(self):
        mu, nu = sp.symbols('mu nu')
        numerator = 'k1(mu)*k2(mu)*k2(1)*k2(2)*p1(1)*p2(2)'
        indices = [mu, nu, 1, 2]
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator, Lorentz_indices=indices)
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [[(0,mu),(1,mu),(1,1),(1,2)]]

        self.assertEqual(len(tensors), 1)
        for i in range(1):
            print(i)
            self.assertEqual(tensors[i], target_tensors[i])

    #@attr('active')
    def test_tensor_numerator_bubble_1L(self):
        numerator = sympify_expression('k(1)*k(2)*k(3)*k(4)')
        contracted_numerator = numerator * sympify_expression('p(1)*p(2)*p(3)*p(4)*const')
        partially_contracted_numerator = numerator * sympify_expression('p(3)*p(4)')
        indices = [1, 2, 3, 4]
        loop_momenta = ['k']
        external_momenta = ['p']
        propagators = ['k**2', '(k - p)**2']

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                           numerator=numerator, Feynman_parameters=['x1','x2'],
                                           Lorentz_indices=indices, dimensionality='D')
        li_contracted = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                                      numerator=contracted_numerator, dimensionality='D',
                                                      Feynman_parameters=['x1','x2'], Lorentz_indices=indices)
        li_partially_contracted = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                                                numerator=partially_contracted_numerator,
                                                                Feynman_parameters=['x1','x2'], Lorentz_indices=indices,
                                                                replacement_rules=[('p*p', 'm**2')], dimensionality='D')
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [[(0,1),(0,2),(0,3),(0,4)]]

        self.assertEqual(tensors, target_tensors)

        numerator = sympify_expression(li.numerator) * li.Gamma_factor
        contracted_numerator = sympify_expression(li_contracted.numerator) * li.Gamma_factor
        partially_contracted_numerator = sympify_expression(li_partially_contracted.numerator) * li.Gamma_factor
        scalar_factor_insertation = { # N_nu = 2, L = 1 in this example
                                          sympify_expression('scalar_factor(0)'): sympify_expression('1/(-2)**(0/2)*(2 - D*1/2 - 2/2)*(2 - D*1/2 - 4/2)*gamma(2 - D*1/2 - 4/2)*F**(0/2)'),
                                          sympify_expression('scalar_factor(2)'): sympify_expression('1/(-2)**(2/2)*(2 - D*1/2 - 4/2)*gamma(2 - D*1/2 - 4/2)*F**(2/2)'),
                                          sympify_expression('scalar_factor(4)'): sympify_expression('1/(-2)**(4/2)*gamma(2 - D*1/2 - 4/2)*F**(4/2)')
                                    }
        target_numerator = sympify_expression('''
                                             scalar_factor(0)*p(1)*x2*p(2)*x2*p(3)*x2*p(4)*x2 +
                                             g(1,2) * scalar_factor(2)*p(3)*x2*p(4)*x2 +
                                             g(1,3) * scalar_factor(2)*p(2)*x2*p(4)*x2 +
                                             g(1,4) * scalar_factor(2)*p(2)*x2*p(3)*x2 +
                                             g(2,3) * scalar_factor(2)*p(1)*x2*p(4)*x2 +
                                             g(2,4) * scalar_factor(2)*p(1)*x2*p(3)*x2 +
                                             g(3,4) * scalar_factor(2)*p(1)*x2*p(2)*x2 +
                                             g(1,2) * g(3,4) * scalar_factor(4) +
                                             g(1,3) * g(2,4) * scalar_factor(4) +
                                             g(1,4) * g(2,3) * scalar_factor(4)
                                      ''').subs(scalar_factor_insertation)
        target_contracted_numerator = sympify_expression('''
                                                        scalar_factor(0)*p(1)**2*x2*p(2)**2*x2*p(3)**2*x2*p(4)**2*x2 +
                                                        p(2)*p(2) * scalar_factor(2)*p(3)*p(3)*x2*p(4)*p(4)*x2 +
                                                        p(3)*p(3) * scalar_factor(2)*p(2)*p(2)*x2*p(4)*p(4)*x2 +
                                                        p(4)*p(4) * scalar_factor(2)*p(2)*p(2)*x2*p(3)*p(3)*x2 +
                                                        p(3)*p(3) * scalar_factor(2)*p(1)*p(1)*x2*p(4)*p(4)*x2 +
                                                        p(4)*p(4) * scalar_factor(2)*p(1)*p(1)*x2*p(3)*p(3)*x2 +
                                                        p(4)*p(4) * scalar_factor(2)*p(1)*p(1)*x2*p(2)*p(2)*x2 +
                                                        p(2)*p(2) * p(4)*p(4) * scalar_factor(4) +
                                                        p(3)*p(3) * p(4)*p(4) * scalar_factor(4) +
                                                        p(4)*p(4) * p(3)*p(3) * scalar_factor(4)
                                                 ''').subs(scalar_factor_insertation) * sympify_expression('const')
        target_partially_contracted_numerator = sympify_expression('''
                                                                  scalar_factor(0)*m**2*x2*m**2*x2*p(1)*x2*p(2)*x2 +
                                                                  m**2 * scalar_factor(2)*p(1)*x2*p(2)*x2 +
                                                                  p(1) * scalar_factor(2)*m**2*x2*p(2)*x2 +
                                                                  p(2) * scalar_factor(2)*m**2*x2*p(1)*x2 +
                                                                  p(1) * scalar_factor(2)*m**2*x2*p(2)*x2 +
                                                                  p(2) * scalar_factor(2)*m**2*x2*p(1)*x2 +
                                                                  g(1,2) * scalar_factor(2)*m**2*x2*m**2*x2 +
                                                                  m**2 * g(1,2) * scalar_factor(4) +
                                                                  p(1) * p(2) * scalar_factor(4) +
                                                                  p(2) * p(1) * scalar_factor(4)
                                                           ''').subs(scalar_factor_insertation)
        self.assertEqual( (numerator - target_numerator).simplify() , 0 )
        self.assertEqual( (contracted_numerator - target_contracted_numerator).simplify() , 0 )
        self.assertEqual( (partially_contracted_numerator - target_partially_contracted_numerator).simplify() , 0 )

    #@attr('active')
    def test_replacement_rules_in_numerator_bubble_1L(self):
        numerator = sympify_expression('k(1)*k(2)*k(3)*k(4)*p(1)*p(2)*p(3)*p(4)')
        indices = [1, 2, 3, 4]
        loop_momenta = ['k']
        external_momenta = ['p']
        propagators = ['k**2', '(k - p)**2']
        replacement_rules = [('p**2', 'm**2')]

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator, Feynman_parameters=['x1','x2'], replacement_rules=replacement_rules, Lorentz_indices=indices)

        target_numerator = sympify_expression('''
                                             scalar_factor(0)*m**2*x2*m**2*x2*m**2*x2*m**2*x2 +
                                             m**2 * scalar_factor(2)*m**2*x2*m**2*x2 +
                                             m**2 * scalar_factor(2)*m**2*x2*m**2*x2 +
                                             m**2 * scalar_factor(2)*m**2*x2*m**2*x2 +
                                             m**2 * scalar_factor(2)*m**2*x2*m**2*x2 +
                                             m**2 * scalar_factor(2)*m**2*x2*m**2*x2 +
                                             m**2 * scalar_factor(2)*m**2*x2*m**2*x2 +
                                             m**2 * m**2 * scalar_factor(4) +
                                             m**2 * m**2 * scalar_factor(4) +
                                             m**2 * m**2 * scalar_factor(4)
                                      ''').subs({ # N_nu = 2, L = 1, D=4-2*eps in this example
                                                      sympify_expression('scalar_factor(0)'): sympify_expression('1/(-2)**(0/2)*(2 - (4-2*eps)*1/2 - 2/2)*(2 - (4-2*eps)*1/2 - 4/2)*gamma(2 - (4-2*eps)*1/2 - 4/2)*F**(0/2)'),
                                                      sympify_expression('scalar_factor(2)'): sympify_expression('1/(-2)**(2/2)*(2 - (4-2*eps)*1/2 - 4/2)*gamma(2 - (4-2*eps)*1/2 - 4/2)*F**(2/2)'),
                                                      sympify_expression('scalar_factor(4)'): sympify_expression('1/(-2)**(4/2)*gamma(2 - (4-2*eps)*1/2 - 4/2)*F**(4/2)')
                                                })

        self.assertEqual( (sympify_expression(li.numerator)*li.Gamma_factor - target_numerator).simplify() , 0 )

    #@attr('active')
    def test_2L_box_with_numerator(self):
        propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2','(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2']
        numerator = '2*k1(mu)*p3(mu) + eps'
        loop_momenta = ['k1','k2']
        external_momenta = ['p1','p2','p3','p4']

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator, Feynman_parameters=['z%i'%i for i in range(1,7+1)], Lorentz_indices=['mu'])

        # comparison with Mathematica implementation of SecDec
        target_U = sympify_expression('''
                                   z4*(z5 + z6 + z7) + z1*(z4 + z5 + z6 + z7) +
                                   z2*(z4 + z5 + z6 + z7) + z3*(z4 + z5 + z6 + z7)
                              ''')
        target_F = sympify_expression('''
                                    z4*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                                    2*p2*p3 + p3*p3)*z6*z7) -
                                    z5*(p1*p1*z6 + 2*p1*p2*z6 + p2*p2*z6 +
                                    p3*p3*z7)) +
                                    z2*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                                    2*p2*p3 + p3*p3)*z6*z7) -
                                    (p1*p1 + 2*p1*p2 + p2*p2)*z3*(z4 + z5 + z6 +
                                    z7) - z4*(p1*p1*z6 + 2*p1*p2*z6 + p2*p2*z6 +
                                    p3*p3*z7) - z5*(p1*p1*z6 + 2*p1*p2*z6 +
                                    p2*p2*z6 + p3*p3*z7)) +
                                    z3*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                                    2*p2*p3 + p3*p3)*z6*z7) -
                                    z5*(p1*p1*z6 + 2*p1*p2*z6 + p2*p2*z6 +
                                    p3*p3*z7) - z4*(p2*p2*z5 + 2*p1*p3*z7 +
                                    p2*p2*z7 + 2*p2*p3*z7 + p3*p3*z7 +
                                    p1*p1*(z5 + z7) + 2*p1*p2*(z5 + z7))) +
                                    z1*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                                    2*p2*p3 + p3*p3)*z6*z7) -
                                    p2*p2*z2*(z4 + z5 + z6 + z7) -
                                    p1*p1*z3*(z4 + z5 + z6 + z7) -
                                    z5*(p1*p1*z6 + 2*p1*p2*z6 + p2*p2*z6 +
                                    p3*p3*z7) - z4*(p1*p1*z6 + (2*p2*p3 + p3*p3)*
                                    z7 + p2*p2*(z5 + z7)))
                              ''')
        target_numerator = sympify_expression('''
                                               2*(-(p3(mu)* p3(mu)*z4*z7)
                                             - p2(mu)* p3(mu)*(z4*(z5 + z7)
                                             + z2*(z4 + z5 + z6 + z7))
                                             + p1(mu)* p3(mu)*(z4*z6 + z3*(z4 + z5 + z6 + z7)))
                                             + eps * U
                                      ''') * sympify_expression('gamma(7 - (4-2*eps)*2/2)')
        # SecDec puts the ``(-1)**N_nu`` into the prefactor while `LoopIntegral` puts it into the `Gamma_factor`.
        target_numerator = - target_numerator

        self.assertEqual( (sympify_expression(li.U) - target_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(li.F) - target_F  ).simplify() , 0 )
        self.assertEqual( (sympify_expression(li.numerator) * li.Gamma_factor - target_numerator).simplify() , 0 )

        replacement_rules = np.array([('p1**2', 0),
                                      ('p2**2', 0),
                                      ('p3**2', 0),
                                      ('p4**2', 0),
                                      ('p1*p2', 's/2'),
                                      ('p2*p3', 't/2'),
                                      ('p1*p3', '-s/2 - t/2')])

        li_with_replacement_rules = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                                                  numerator=numerator, Feynman_parameters=['z%i'%i for i in range(1,7+1)],
                                                                  replacement_rules=replacement_rules, Lorentz_indices=['mu'])

        target_numerator_with_replacements = sympify_expression('''
                                                             2*(
                                                                    - (t/2)*(z4*(z5 + z7)
                                                                    + z2*(z4 + z5 + z6 + z7))
                                                                    + (-s/2-t/2)*(z4*z6 + z3*(z4 + z5 + z6 + z7))
                                                               ) + eps * U
                                                        ''') * sympify_expression('gamma(7 - (4-2*eps)*2/2)')
        # SecDec puts the ``(-1)**N_nu`` into the prefactor while `LoopIntegral` puts it into the `Gamma_factor`.
        target_numerator_with_replacements = - target_numerator_with_replacements

        self.assertEqual( (sympify_expression(li_with_replacement_rules.numerator) * li.Gamma_factor - target_numerator_with_replacements).simplify() , 0 )

    #@attr('active')
    @attr('slow')
    def test_rank3_numerator_2L_double_triangle(self):
        k1,  k2,  p1,  p2,  p3,  p4, s, t, mu, nu, scalar_factor, eps, F = sp.symbols('k1 k2 p1 p2 p3 p4 s t mu nu scalar_factor eps F')
        k1f, k2f, p1f, p2f, p3f, p4f = sp.symbols('k1 k2 p1 p2 p3 p4', cls=sp.Function)
        sp11, sp12, sp13, sp22, sp23, sp33 = sp.symbols('sp11 sp12 sp13 sp22 sp23 sp33')

        z = sympify_expression(['z%i'%i for i in range(6+1)]) # Feynman parameters (only use z[1], z[2], ..., z[6] but not z[0])
        loop_momenta = [k1, k2]
        external_momenta = [p1, p2, p3, p4]
        indices = [mu, nu]
        replacement_rules = [('p4', '- (p1 + p2 + p3)'),
                             ('p1*p1', 'sp11'),
                             ('p2*p2', 'sp22'),
                             ('p3*p3', 'sp33'),
                             ('p1*p2', 'sp12'),
                             ('p2*p3', 'sp23'),
                             ('p1*p3', 'sp13')]

        propagators = [k1**2,(k1+p1)**2,(k1-k2-p3-p4)**2,(k2)**2,(k2+p2)**2,(k2+p3)**2]
        numerator = k1f(mu)*k2f(nu)*k1f(mu)*p2f(nu)*2

        li = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                           numerator=numerator, Feynman_parameters=z[1:],
                                           Lorentz_indices=indices,
                                           replacement_rules=replacement_rules)

        # need to insert the `scalar_factor` (and the `F` therein) when comparing with SecDec
        # `scalar_factor` = `1/(-2)**(r/2)*Gamma(N_nu - dim*L/2 - r/2)*F**(r/2)
        N_nu = sympify_expression(len(propagators))
        dim = sympify_expression(4 - 2*eps)
        L = sympify_expression(2)
        numerator = sympify_expression(li.numerator) * li.Gamma_factor
        numerator = numerator.subs(F, sympify_expression(li.F))

        # SecDec divides by a factor of ``Gamma(N_nu - dim*L/2)`` and uses ``Gamma(N_nu - dim*L/2 - 1) = Gamma(N_nu - dim*L/2) / (N_nu - dim*L/2 - 1)``
        numerator /= sp.gamma(N_nu - dim*L/2)
        numerator  = numerator.subs(sp.gamma(N_nu - dim*L/2 - 1), sp.gamma(N_nu - dim*L/2) / (N_nu - dim*L/2 - 1))


        # Results of the Mathematica implementation
        target_U = z[3]*(z[4] + z[5] + z[6]) + z[1]*(z[3] + z[4] + z[5] + z[6]) + \
                   z[2]*(z[3] + z[4] + z[5] + z[6])
        target_F = z[3]*(-((-2*sp23*z[5] + sp33*(z[4] + z[5]))*z[6]) - \
                   sp22*z[5]*(z[4] + z[6])) + \
                   z[1]*(-((sp22 - 2*sp23 + sp33)*z[5]*z[6]) - \
                   sp11*z[2]*(z[3] + z[4] + z[5] + z[6]) - \
                   z[4]*(sp22*z[5] + sp33*z[6]) - \
                   z[3]*(sp22*z[4] + 4*sp22*z[5] + 2*sp13*z[6] + \
                   sp22*z[6] + 2*sp23*z[6] + sp33*z[6] + \
                   sp11*(z[4] + z[5] + z[6]) + 2*sp12*\
                   (z[4] + 2*z[5] + z[6]))) + \
                   z[2]*(-((2*sp23*(z[3] - z[5]) + sp33*(z[3] + z[4] + z[5]))*\
                   z[6]) - sp22*(z[5]*(z[4] + z[6]) + z[3]*(z[4] + 4*z[5] + z[6])))
        target_numerator = 2*(-((4 - 2*eps)*(z[3] + z[4] + z[5] + z[6])*(sp12*z[1]*z[3] + \
                           sp22*(z[1]*(z[3] - z[5]) + z[2]*(z[3] - z[5]) - z[3]*z[5]) - \
                           sp23*(z[1] + z[2] + z[3])*z[6])*\
                           (z[3]*(-((-2*sp23*z[5] + sp33*(z[4] + z[5]))*z[6]) - \
                           sp22*z[5]*(z[4] + z[6])) + \
                           z[1]*(-((sp22 - 2*sp23 + sp33)*z[5]*z[6]) - \
                           sp11*z[2]*(z[3] + z[4] + z[5] + z[6]) - \
                           z[4]*(sp22*z[5] + sp33*z[6]) - \
                           z[3]*(sp22*z[4] + 4*sp22*z[5] + 2*sp13*z[6] + \
                           sp22*z[6] + 2*sp23*z[6] + sp33*z[6] + \
                           sp11*(z[4] + z[5] + z[6]) + 2*sp12*(z[4] + 2*z[5] + \
                           z[6]))) + z[2]*(-((2*sp23*(z[3] - z[5]) + \
                           sp33*(z[3] + z[4] + z[5]))*z[6]) - \
                           sp22*(z[5]*(z[4] + z[6]) + z[3]*(z[4] + 4*z[5] + z[6])))))/ \
                           (2*(1 + 2*eps)) - \
                           (z[3]*(-(sp12*(z[3]*(z[4] + z[5] + z[6]) + \
                           z[2]*(z[3] + z[4] + z[5] + z[6]))) - \
                           z[3]*(sp23*z[6] + sp22*(z[4] + 2*z[5] + z[6])))*\
                           (z[3]*(-((-2*sp23*z[5] + sp33*(z[4] + z[5]))*z[6]) - \
                           sp22*z[5]*(z[4] + z[6])) + \
                           z[1]*(-((sp22 - 2*sp23 + sp33)*z[5]*z[6]) - \
                           sp11*z[2]*(z[3] + z[4] + z[5] + z[6]) - \
                           z[4]*(sp22*z[5] + sp33*z[6]) - \
                           z[3]*(sp22*z[4] + 4*sp22*z[5] + 2*sp13*z[6] + \
                           sp22*z[6] + 2*sp23*z[6] + sp33*z[6] + \
                           sp11*(z[4] + z[5] + z[6]) + 2*sp12*(z[4] + 2*z[5] + \
                           z[6]))) + z[2]*(-((2*sp23*(z[3] - z[5]) + \
                           sp33*(z[3] + z[4] + z[5]))*z[6]) - \
                           sp22*(z[5]*(z[4] + z[6]) + z[3]*(z[4] + 4*z[5] + z[6])))))/ \
                           (1 + 2*eps) + (sp12*z[1]*z[3] + \
                           sp22*(z[1]*(z[3] - z[5]) + z[2]*(z[3] - z[5]) - z[3]*z[5]) - \
                           sp23*(z[1] + z[2] + z[3])*z[6])*\
                           (sp11*(z[3]*(z[4] + z[5] + z[6]) + z[2]*(z[3] + z[4] + z[5] + z[6]))**\
                           2 + z[3]*(sp22*z[3]*(z[4] + 2*z[5] + z[6])**2 + \
                           2*sp12*(z[4] + 2*z[5] + z[6])*(z[3]*(z[4] + z[5] + z[6]) + \
                           z[2]*(z[3] + z[4] + z[5] + z[6])) + \
                           z[6]*(2*sp13*(z[3]*(z[4] + z[5] + z[6]) + \
                           z[2]*(z[3] + z[4] + z[5] + z[6])) + \
                           z[3]*(sp33*z[6] + 2*sp23*(z[4] + 2*z[5] + z[6]))))))

        self.assertEqual( (sympify_expression(li.U) - target_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(li.F) - target_F).simplify() , 0 )
        self.assertEqual( (sympify_expression(numerator) - target_numerator).simplify() , 0 )


        rank = 3
        target_exponent_U = N_nu - dim/2 * (L+1) - rank
        target_exponent_F = -(N_nu - dim/2 * L)
        target_exponentiated_F = target_F ** target_exponent_F
        target_exponentiated_U = target_U ** target_exponent_U
        self.assertEqual( (sympify_expression(li.exponentiated_U) - target_exponentiated_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(li.exponentiated_F) - target_exponentiated_F).simplify() , 0 )

    #@attr('active')
    @attr('slow')
    def test_rank2_numerator_2L_box(self):
        k1,  k2,  p1,  p2,  p3,  p4, s, t, mu, scalar_factor, D, eps = sp.symbols('k1 k2 p1 p2 p3 p4 s t mu scalar_factor D eps')
        k1f, k2f, p1f, p2f, p3f, p4f = sp.symbols('k1 k2 p1 p2 p3 p4', cls=sp.Function)

        z = sympify_expression(['z%i'%i for i in range(7+1)]) # Feynman parameters (only use z[1], z[2], ..., z[7] but not z[0])
        loop_momenta = [k1, k2]
        external_momenta = [p1, p2, p3, p4]
        indices = [mu]

        propagators = [k1**2,(k1+p2)**2,(k1-p1)**2,(k1-k2)**2,(k2+p2)**2,(k2-p1)**2,(k2+p2+p3)**2]
        numerator1 = 2*k1f(mu)*k2f(mu)
        numerator2 = 2*k1f(mu)*k1f(mu)
        replacement_rules = [(p1**2, 0),
                             (p2**2, 0),
                             (p3**2, 0),
                             (p4**2, 0),
                             (p1*p2, s/2),
                             (p2*p3, t/2),
                             (p1*p3, -s/2 - t/2)]

        li1 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                            numerator=numerator1, Feynman_parameters=z[1:],
                                            Lorentz_indices=indices, dimensionality=D,
                                            replacement_rules=replacement_rules)
        li2 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta,
                                            numerator=numerator2, Feynman_parameters=z[1:],
                                            Lorentz_indices=indices)
        Feynman_parametrized_numerator1 = sympify_expression(li1.numerator)
        Feynman_parametrized_numerator2 = sympify_expression(li2.numerator)

        # comparison with Mathematica implementation of SecDec
        target_U = z[4]*(z[5] + z[6] + z[7]) + z[1]*(z[4] + z[5] + z[6] + z[7]) + \
                   z[2]*(z[4] + z[5] + z[6] + z[7]) + z[3]*(z[4] + z[5] + z[6] + z[7])
        target_F1 = -(t*z[1]*z[4]*z[7]) - \
                    s*(z[5]*((z[1] + z[4])*z[6] + z[3]*(z[4] + z[6])) + \
                    z[2]*((z[4] + z[5])*z[6] + z[3]*(z[4] + z[5] + z[6] + z[7])))
        target_F2 = z[4]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 + \
                    2*p2*p3 + p3*p3)*z[6]*z[7]) - \
                    z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] + \
                    p3*p3*z[7])) + \
                    z[2]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 + \
                    2*p2*p3 + p3*p3)*z[6]*z[7]) - \
                    (p1*p1 + 2*p1*p2 + p2*p2)*z[3]*(z[4] + z[5] + z[6] + \
                    z[7]) - z[4]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] + \
                    p3*p3*z[7]) - z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + \
                    p2*p2*z[6] + p3*p3*z[7])) + \
                    z[3]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 + \
                    2*p2*p3 + p3*p3)*z[6]*z[7]) - \
                    z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] + \
                    p3*p3*z[7]) - z[4]*(p2*p2*z[5] + 2*p1*p3*z[7] + \
                    p2*p2*z[7] + 2*p2*p3*z[7] + p3*p3*z[7] + \
                    p1*p1*(z[5] + z[7]) + 2*p1*p2*(z[5] + z[7]))) + \
                    z[1]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 + \
                    2*p2*p3 + p3*p3)*z[6]*z[7]) - \
                    p2*p2*z[2]*(z[4] + z[5] + z[6] + z[7]) - \
                    p1*p1*z[3]*(z[4] + z[5] + z[6] + z[7]) - \
                    z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] + \
                    p3*p3*z[7]) - z[4]*(p1*p1*z[6] + (2*p2*p3 + p3*p3)* \
                    z[7] + p2*p2*(z[5] + z[7])))
        target_Fs = [target_F1, target_F2]
        target_numerator1 = (-(s*((-6 + D)*z[2]**2*z[6]*(z[4] + z[5] + z[6] + z[7]) + \
                            z[5]*(3*(-4 + D)*z[4]**2*z[6] + (-6 + D)*z[3]**2* \
                            (z[4] + z[5] + z[6] + z[7]) + z[3]*z[4]* \
                            (3*(-4 + D)*z[4] + (-6 + D)*z[5] - \
                            18*z[6] + 4*D*z[6] - 6*z[7] + D*z[7]) + \
                            z[1]*(3*(-4 + D)*z[4]*z[6] + (-6 + D)*z[3]* \
                            (z[4] + z[5] + z[6] + z[7]))) + \
                            z[2]*(z[3]*(3*(-4 + D)*z[4] + (-6 + D)* \
                            (z[5] + z[6]))*(z[4] + z[5] + z[6] + z[7]) + \
                            z[6]*((-6 + D)*z[1]*(z[4] + z[5] + z[6] + z[7]) + \
                            z[4]*(3*(-4 + D)*z[4] + 2*(-9 + 2*D)* \
                            z[5] + (-6 + D)*(z[6] + z[7])))))) + \
                            t*z[7]*((-6 + D)*(z[2] + z[3] + 2*z[4])* \
                            (z[4]*(z[5] + z[6] + z[7]) + z[2]*(z[4] + z[5] + z[6] + z[7]) + \
                            z[3]*(z[4] + z[5] + z[6] + z[7])) + \
                            z[1]*((-6 + D)*z[2]*(z[4] + z[5] + z[6] + z[7]) + \
                            (-6 + D)*z[3]*(z[4] + z[5] + z[6] + z[7]) - \
                            z[4]*(12*(z[5] + z[6] + z[7]) + D* \
                            (z[4] - 2*(z[5] + z[6] + z[7]))))))/(-6 + D)
        target_numerator2 = (2*(p2*p2*z[2]**2*z[4]**2 - 2*p1*p2*z[2]*z[3]*z[4]**2 +
                            p1*p1*z[3]**2*z[4]**2 + 2*p2*p2*z[2]**2*z[4]*z[5] -
                            4*p1*p2*z[2]*z[3]*z[4]*z[5] + 2*p1*p1*z[3]**2*z[4]*z[5] +
                            2*p2*p2*z[2]*z[4]**2*z[5] - 2*p1*p2*z[3]*z[4]**2*z[5] +
                            p2*p2*z[2]**2*z[5]**2 - 2*p1*p2*z[2]*z[3]*z[5]**2 +
                            p1*p1*z[3]**2*z[5]**2 + 2*p2*p2*z[2]*z[4]*z[5]**2 -
                            2*p1*p2*z[3]*z[4]*z[5]**2 + p2*p2*z[4]**2*z[5]**2 +
                            2*p2*p2*z[2]**2*z[4]*z[6] - 4*p1*p2*z[2]*z[3]*z[4]*z[6] +
                            2*p1*p1*z[3]**2*z[4]*z[6] - 2*p1*p2*z[2]*z[4]**2*z[6] +
                            2*p1*p1*z[3]*z[4]**2*z[6] + 2*p2*p2*z[2]**2*z[5]*z[6] -
                            4*p1*p2*z[2]*z[3]*z[5]*z[6] + 2*p1*p1*z[3]**2*z[5]*z[6] -
                            2*p1*p2*z[2]*z[4]*z[5]*z[6] + 2*p2*p2*z[2]*z[4]*z[5]*z[6] +
                            2*p1*p1*z[3]*z[4]*z[5]*z[6] - 2*p1*p2*z[3]*z[4]*z[5]*z[6] -
                            2*p1*p2*z[4]**2*z[5]*z[6] + p2*p2*z[2]**2*z[6]**2 -
                            2*p1*p2*z[2]*z[3]*z[6]**2 + p1*p1*z[3]**2*z[6]**2 -
                            2*p1*p2*z[2]*z[4]*z[6]**2 + 2*p1*p1*z[3]*z[4]*z[6]**2 +
                            p1*p1*z[4]**2*z[6]**2 + 2*p2*p2*z[2]**2*z[4]*z[7] -
                            4*p1*p2*z[2]*z[3]*z[4]*z[7] + 2*p1*p1*z[3]**2*z[4]*z[7] +
                            2*p2*p2*z[2]*z[4]**2*z[7] + 2*p2*p3*z[2]*z[4]**2*z[7] -
                            2*p1*p2*z[3]*z[4]**2*z[7] - 2*p1*p3*z[3]*z[4]**2*z[7] +
                            2*p2*p2*z[2]**2*z[5]*z[7] - 4*p1*p2*z[2]*z[3]*z[5]*z[7] +
                            2*p1*p1*z[3]**2*z[5]*z[7] + 4*p2*p2*z[2]*z[4]*z[5]*z[7] +
                            2*p2*p3*z[2]*z[4]*z[5]*z[7] - 4*p1*p2*z[3]*z[4]*z[5]*z[7] -
                            2*p1*p3*z[3]*z[4]*z[5]*z[7] + 2*p2*p2*z[4]**2*z[5]*z[7] +
                            2*p2*p3*z[4]**2*z[5]*z[7] + 2*p2*p2*z[2]**2*z[6]*z[7] -
                            4*p1*p2*z[2]*z[3]*z[6]*z[7] + 2*p1*p1*z[3]**2*z[6]*z[7] -
                            2*p1*p2*z[2]*z[4]*z[6]*z[7] + 2*p2*p2*z[2]*z[4]*z[6]*z[7] +
                            2*p2*p3*z[2]*z[4]*z[6]*z[7] + 2*p1*p1*z[3]*z[4]*z[6]*z[7] -
                            2*p1*p2*z[3]*z[4]*z[6]*z[7] - 2*p1*p3*z[3]*z[4]*z[6]*z[7] -
                            2*p1*p2*z[4]**2*z[6]*z[7] - 2*p1*p3*z[4]**2*z[6]*z[7] +
                            p2*p2*z[2]**2*z[7]**2 - 2*p1*p2*z[2]*z[3]*z[7]**2 +
                            p1*p1*z[3]**2*z[7]**2 + 2*p2*p2*z[2]*z[4]*z[7]**2 +
                            2*p2*p3*z[2]*z[4]*z[7]**2 - 2*p1*p2*z[3]*z[4]*z[7]**2 -
                            2*p1*p3*z[3]*z[4]*z[7]**2 + p2*p2*z[4]**2*z[7]**2 +
                            2*p2*p3*z[4]**2*z[7]**2 + p3*p3*z[4]**2*z[7]**2 -
                            ((4 - 2*eps)*(z[4] + z[5] + z[6] + z[7])*
                            (z[4]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                            2*p2*p3 + p3*p3)*z[6]*z[7]) -
                            z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] +
                            p3*p3*z[7])) +
                            z[2]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                            2*p2*p3 + p3*p3)*z[6]*z[7]) -
                            (p1*p1 + 2*p1*p2 + p2*p2)*z[3]*(z[4] + z[5] + z[6] +
                            z[7]) - z[4]*(p1*p1*z[6] + 2*p1*p2*z[6] +
                            p2*p2*z[6] + p3*p3*z[7]) -
                            z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] +
                            p3*p3*z[7])) +
                            z[3]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                            2*p2*p3 + p3*p3)*z[6]*z[7]) -
                            z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] + p2*p2*z[6] +
                            p3*p3*z[7]) - z[4]*(p2*p2*z[5] + 2*p1*p3*z[7] +
                            p2*p2*z[7] + 2*p2*p3*z[7] + p3*p3*z[7] +
                            p1*p1*(z[5] + z[7]) + 2*p1*p2*(z[5] + z[7]))) +
                            z[1]*(-((p1*p1 + 2*p1*p2 + 2*p1*p3 + p2*p2 +
                            2*p2*p3 + p3*p3)*z[6]*z[7]) - p2*p2*z[2]*
                            (z[4] + z[5] + z[6] + z[7]) - p1*p1*z[3]*(z[4] + z[5] + z[6] +
                            z[7]) - z[5]*(p1*p1*z[6] + 2*p1*p2*z[6] +
                            p2*p2*z[6] + p3*p3*z[7]) -
                            z[4]*(p1*p1*z[6] + (2*p2*p3 + p3*p3)*z[7] +
                            p2*p2*(z[5] + z[7])))))/(4*(1 + eps))))
        target_numerators = [target_numerator1, target_numerator2]

        self.assertEqual( (sympify_expression(li1.U) - target_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(li2.U) - target_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(li1.F) - target_F1).simplify() , 0 )
        self.assertEqual( (sympify_expression(li2.F) - target_F2).simplify() , 0 )

        for i, Feynman_parametrized_numerator in enumerate([Feynman_parametrized_numerator1, Feynman_parametrized_numerator2]):
            print(i)

            # need to insert the `F`s when comparing with SecDec
            Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs('F', target_Fs[i])

            # SecDec divides by a factor of ``Gamma(7 - D)`` while pySecDec divides by ``Gamma(6 - D)``
            Feynman_parametrized_numerator *= li1.Gamma_factor / sympify_expression('gamma(7 - D)')
            Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs('gamma(7 - D)', '(6 - D)*gamma(6 - D)')

            # SecDec puts the ``(-1)**N_nu`` into the prefactor while `LoopIntegral` puts it into the `Gamma_factor`.
            Feynman_parametrized_numerator = - Feynman_parametrized_numerator

            # only one scalar product per term --> can safely remove dummy index
            Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs(p1f(mu), p1).subs(p2f(mu), p2).subs(p3f(mu), p3).subs(p4f(mu), p4)

            if i == 1:
                Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs(D, 4-2*eps)

            self.assertEqual( (sympify_expression(Feynman_parametrized_numerator) - target_numerators[i]).simplify() , 0 )

    #@attr('active')
    def test_numerator_without_external_momenta(self):
        double_tadpole = LoopIntegralFromPropagators(
            loop_momenta = ['k1','k2'],
            Lorentz_indices = ['mu'],
            propagators = ['(k1)**2-msq','k2**2'],
            numerator = 'k1(mu)*k1(mu)',
        )

        D,eps = sp.symbols('D eps')

        target_F = sympify_expression('msq * x0**2*x1')
        target_exponent_F = sympify_expression('-(Nnu-L*D/2)').subs(D,4-2*eps).subs('R',2).subs('L',2).subs('Nnu',2)
        target_U = sympify_expression('x0*x1')
        target_exponent_U = sympify_expression('Nnu-(L+1)*D/2-R').subs(D,4-2*eps).subs('R',2).subs('L',2).subs('Nnu',2)
        target_numerator = sympify_expression('-D/2 * x1 * F').subs(D,4-2*eps)

        self.assertEqual( (sympify_expression(double_tadpole.F) - target_F).simplify() , 0 )
        self.assertEqual( (sympify_expression(double_tadpole.exponent_F) - target_exponent_F).simplify() , 0 )
        self.assertEqual( (sympify_expression(double_tadpole.U) - target_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(double_tadpole.exponent_U) - target_exponent_U).simplify() , 0 )
        self.assertEqual( (sympify_expression(double_tadpole.numerator) - target_numerator).simplify() , 0 )

    #@attr('active')
    def test_error_if_index_too_often_loop_momenta(self):
        li = LoopIntegralFromPropagators(
            Lorentz_indices = ('mu', 1, 2),
            numerator = 'k1(1)*k1(1)*k1(1) + k2(mu)*k2(2)*p1(mu)*p2(2)',
            loop_momenta = ['k1', 'k2'],
            external_momenta = ['p1', 'p2'],
            propagators = ['(k1-p1)^2','(k2-p1-p2)^2']
        )
        self.assertRaisesRegexp(
            AssertionError,
            str(li.numerator_input.args[0]).replace('(',r'\(').replace(')','\)').replace('*',r'\*') + \
                '.*Lorentz_indices.*1.*more than.*(twice|two times)',
            lambda: li.numerator
        )

    #@attr('active')
    def test_error_if_index_too_often_external_momenta(self):
        li = LoopIntegralFromPropagators(
            Lorentz_indices = ('mu', 1, 2),
            numerator = 'p1(1)*p1(1)*p2(1) + k2(mu)*k2(2)*p1(mu)*p2(2)',
            loop_momenta = ['k1', 'k2'],
            external_momenta = ['p1', 'p2'],
            propagators = ['(k1-p1)^2','(k2-p1-p2)^2']
        )
        self.assertRaisesRegexp(
            AssertionError,
            str(li.numerator_input.args[0]).replace('(',r'\(').replace(')','\)').replace('*',r'\*') + \
                '.*Lorentz_indices.*1.*more than.*(twice|two times)',
            lambda: li.numerator
        )

    #@attr('active')
    def test_error_if_index_too_often_mixed(self):
        li = LoopIntegralFromPropagators(
            Lorentz_indices = ('mu', 1, 2),
            numerator = 'k1(mu)*k2(mu)*k2(1)*k2(2)*p1(1)*p2(2)*p1(mu)',
            loop_momenta = ['k1', 'k2'],
            external_momenta = ['p1', 'p2'],
            propagators = ['(k1-p1)^2','(k2-p1-p2)^2']
        )
        self.assertRaisesRegexp(
            AssertionError,
            str(li.numerator_input).replace('(',r'\(').replace(')','\)').replace('*',r'\*') + \
                '.*Lorentz_indices.*mu.*more than.*(twice|two times)',
            lambda: li.numerator
        )


def uf_from_graph_generic(test_case, int_lines, ext_lines, result_L, result_u, result_f, rules=[]):
    loop_integral = LoopIntegralFromGraph(int_lines, ext_lines, Feynman_parameters='Feynman',
                                            replacement_rules=rules)

    test_case.assertEqual(loop_integral.L,result_L)

    u,f = loop_integral.U, loop_integral.F

    sympy_u = sympify_expression(str(u))
    sympy_f = sympify_expression(str(f))

    result_u = sympify_expression(str(result_u).replace('x','Feynman'))
    result_f = sympify_expression(str(result_f).replace('x','Feynman'))

    zerou = (sympy_u - result_u).expand()
    zerof = (sympy_f - result_f).expand()

    test_case.assertEqual(zerou,0)
    test_case.assertEqual(zerof,0)

    expo_u, expo_f = loop_integral.exponent_U, loop_integral.exponent_F
    test_case.assertEqual(expo_u, len(int_lines) - sympify_expression('2-eps') * (1 + loop_integral.L))
    test_case.assertEqual(expo_f, -(len(int_lines) - sympify_expression('2-eps') * loop_integral.L))

#@attr('active')
class TestUF_FromGraph(unittest.TestCase):
    def test_tadpole_1l(self):
        uf_from_graph_generic(self,
                              int_lines = [['m',[1,1]]],
                              ext_lines = [],
                              result_L = 1,
                              result_u = "x0",
                              result_f = "m**2*x0**2")

    #@attr('active')
    def test_bubble_1l(self):
        uf_from_graph_generic(self,
                              int_lines = [['m',(1,2)], [0,(1,2)]],
                              ext_lines = [('p',1),('p',2)],
                              result_L = 1,
                              result_u = "x0 + x1",
                              result_f = "(m**2 - p**2)*x0*x1 + (m**2)*x0**2")

    def test_triangle_1l(self):
        uf_from_graph_generic(self,
                              int_lines = [[0,['vertex1','vertex2']], [0,['vertex2','vertex3']],
                                           [0,['vertex3','vertex1']]],
                              ext_lines = [['p1','vertex1'], ['p2','vertex2'], ['-p1-p2','vertex3']],
                              result_L = 1,
                              result_u = "x0 + x1 + x2",
                              result_f = "-p2sqr*x0*x1 - p1sqr*x0*x2 - p12sqr*x1*x2",
                              rules=[('p1*p1','p1sqr'),
                                     ('p2*p2','p2sqr'),
                                     ('p1*p2','(p12sqr - p1sqr - p2sqr)/2')]
                              )

    def test_bubble_3l(self):
        uf_from_graph_generic(self,
                              int_lines = [[0,[1,2]], [0,[1,2]], [0,[1,2]], [0,[1,2]]],
                              ext_lines = [['p1',1], ['p2',2]],
                              rules = [('p2','-p1')],
                              result_L = 3,
                              result_u = "x0*x1*x2 + x0*x1*x3 + x0*x2*x3 + x1*x2*x3",
                              result_f = "-p1**2*x0*x1*x2*x3")

    # from SecDec -> loop/demos/3_nonplanarbox_2L
    def test_nonplanarbox_2l(self):
        uf_from_graph_generic(self,
                              int_lines=[['m',[1,5]],['m',[2,6]],['M',[1,2]],['M',[3,5]],['m',[3,6]],
                                         ['m',[4,6]],['M',[4,5]]],
                              ext_lines=[['p1',1],['p2',2],['p3',3],['p4',4]],
                              rules=[('p1*p1','msq'),
                                     ('p2*p2','msq'),
                                     ('p3*p3','msq'),
                                     ('p4*p4','msq'),
                                     ('p3*p2','t/2-msq'),
                                     ('p1*p3','-t/2-s/2+msq'),
                                     ('p1*p2','s/2-msq'),
                                     ('p1*p4','t/2-msq'),
                                     ('p2*p4','-t/2-s/2+msq'),
                                     ('p3*p4','s/2-msq'),
                                     ('m**2','msq'),
                                     ('M**2','Msq')],
                              result_L = 2,
                              result_f = """msq*x0**2*x3 + 2*msq*x0*x1*x3 - s*x0*x1*x3 + msq*x1**2*x3
                              + Msq*x0*x2*x3 + Msq*x1*x2*x3 + Msq*x2**2*x3 + Msq*x0*x3**2 + Msq*x1*x3**2
                              + Msq*x2*x3**2 + msq*x0**2*x4 + 2*msq*x0*x1*x4 - s*x0*x1*x4 + msq*x1**2*x4
                              + Msq*x0*x2*x4 + Msq*x1*x2*x4 + Msq*x2**2*x4 + Msq*x0*x3*x4 + Msq*x1*x3*x4
                              + Msq*x2*x3*x4 + msq*x0*x4**2 + msq*x1*x4**2 + msq*x2*x4**2 + msq*x0**2*x5
                              + 2*msq*x0*x1*x5 - s*x0*x1*x5 + msq*x1**2*x5 + Msq*x0*x2*x5 + Msq*x1*x2*x5
                              + Msq*x2**2*x5 + msq*x0*x3*x5 + Msq*x0*x3*x5 + msq*x1*x3*x5 + Msq*x1*x3*x5
                              + msq*x2*x3*x5 + 2*Msq*x2*x3*x5 - t*x2*x3*x5 + Msq*x3**2*x5 + 3*msq*x0*x4*x5
                              - s*x0*x4*x5 + 3*msq*x1*x4*x5 + msq*x2*x4*x5 + Msq*x2*x4*x5 + Msq*x3*x4*x5
                              + msq*x4**2*x5 + msq*x0*x5**2 + msq*x1*x5**2 + msq*x2*x5**2 + msq*x3*x5**2
                              + msq*x4*x5**2 + msq*x0**2*x6 + 2*msq*x0*x1*x6 - s*x0*x1*x6 + msq*x1**2*x6
                              + Msq*x0*x2*x6 + Msq*x1*x2*x6 + Msq*x2**2*x6 + msq*x0*x3*x6 + 2*Msq*x0*x3*x6
                              + msq*x1*x3*x6 + 2*Msq*x1*x3*x6 - s*x1*x3*x6 - msq*x2*x3*x6 + 3*Msq*x2*x3*x6
                              + Msq*x3**2*x6 + msq*x0*x4*x6 + Msq*x0*x4*x6 + msq*x1*x4*x6 + Msq*x1*x4*x6
                              - 3*msq*x2*x4*x6 + 2*Msq*x2*x4*x6 + s*x2*x4*x6 + t*x2*x4*x6 + Msq*x3*x4*x6
                              + msq*x4**2*x6 + Msq*x0*x5*x6 + Msq*x1*x5*x6 + Msq*x2*x5*x6 + Msq*x3*x5*x6
                              + Msq*x4*x5*x6 + Msq*x0*x6**2 + Msq*x1*x6**2 + Msq*x2*x6**2 + Msq*x3*x6**2
                              + Msq*x4*x6**2""",
                              result_u = """x0*x3 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x2*x4 + x0*x5 + x1*x5
                              + x2*x5 + x3*x5 + x4*x5 + x0*x6 + x1*x6 + x2*x6 + x3*x6 + x4*x6"""
                              )

    # from SecDec -> loop/demos/5_pentagon_2L
    def test_pentagon_2l(self):
        uf_from_graph_generic(self,
                              int_lines=[[0,[1,2]],[0,[2,6]],[0,[6,3]],
                                         [0,[3,4]],[0,[4,5]],[0,[5,7]],
                                         [0,[7,1]],[0,[7,6]]],
                              ext_lines = [['p1',1],['p2',2],['p3',3],['p4',4],['p5',5]],
                              rules = [('p1*p1', '0'),
                                       ('p2*p2', '0'),
                                       ('p3*p3', '0'),
                                       ('p4*p4', '0'),
                                       ('p5*p5', '0'),
                                       ('p1*p2', 's12/2'),
                                       ('p1*p3', '(s45-s12-s23)/2'),
                                       ('p1*p4', '(s23-s51-s45)/2'),
                                       ('p1*p5', 's51/2'),
                                       ('p2*p3', 's23/2'),
                                       ('p2*p4', '(-s23-s34+s51)/2'),
                                       ('p2*p5', '(s34-s12-s51)/2'),
                                       ('p3*p4', 's34/2'),
                                       ('p3*p5', '(s12-s34-s45)/2'),
                                       ('p4*p5', 's45/2')],
                              result_L = 2,
                              result_u = """x0*x2 + x1*x2 + x0*x3 + x1*x3 + x0*x4 + x1*x4 + x0*x5 + x1*x5
                              + x2*x6 + x3*x6 + x4*x6 + x5*x6 + x0*x7 + x1*x7 + x2*x7 + x3*x7 + x4*x7
                              + x5*x7 + x6*x7""",
                              result_f = """-s34*x0*x2*x4 - s34*x1*x2*x4 - s12*x0*x2*x5 - s12*x1*x2*x5
                              - s45*x0*x3*x5 - s45*x1*x3*x5 - s12*x1*x2*x6 - s12*x1*x3*x6 - s12*x1*x4*x6
                              - s34*x2*x4*x6 - s12*x1*x5*x6 - s12*x2*x5*x6 - s45*x3*x5*x6 - s23*x0*x3*x7
                              - s51*x0*x4*x7 - s34*x1*x4*x7 - s34*x2*x4*x7 - s12*x1*x5*x7 - s12*x2*x5*x7
                              - s45*x3*x5*x7 - s12*x1*x6*x7 - s12*x2*x6*x7 - s45*x3*x6*x7"""
                              )

    def test_error_messages(self):
        self.assertRaisesRegexp(AssertionError,
                                "To define a loop integral please input a graph with at least one internal line.",
                                LoopIntegralFromGraph, internal_lines = [], external_lines = [['p1',1]])
        self.assertRaisesRegexp(AssertionError,
                                "To define a loop integral please input a graph with at least one closed loop.",
                                LoopIntegralFromGraph, internal_lines = [['m',[1,2]]],
                                external_lines = [['p1',1],['p2',2]])
        self.assertRaisesRegexp(AssertionError,
                                '.*(I|i)nternal.*lines.*form', LoopIntegralFromGraph,
                                internal_lines = [['m',1,1]], external_lines = [])
        self.assertRaisesRegexp(AssertionError,
                                '.*mass.*symbol', LoopIntegralFromGraph,
                                internal_lines = [['m1+m2',[1,1]]], external_lines = [])
        self.assertRaisesRegexp(AssertionError,
                                '.*vertices.*symbol', LoopIntegralFromGraph,
                                internal_lines = [['m',['sin(x)',1]]], external_lines = [])
        self.assertRaisesRegexp(AssertionError,
                                '.*vertices.*symbol', LoopIntegralFromGraph,
                                internal_lines = [['m',[1,1]]], external_lines = [['p1','cos(x)']])
        self.assertRaisesRegexp(AssertionError,
                                '.*propagator.*powers.*vanishing.*regulator', LoopIntegralFromGraph,
                                internal_lines = [['m',[1,1]]], external_lines = [['p1',1]], powerlist=['a+eps'])
        self.assertRaisesRegexp(AssertionError,
                                ".external.*line.*linear.*combination",
                                LoopIntegralFromGraph, internal_lines = [['m',[1,2]], ['m',[1,2]]],
                                external_lines = [['p1',1],['p1**2',2]])

def compare_two_loop_integrals(testcase, li1, li2):
    result_U_1 = sympify_expression(li1.U)
    result_F_1 = sympify_expression(li1.F)
    result_Nu_1 = sympify_expression(li1.numerator).subs('U',result_U_1).subs('F',result_F_1)\
                      *sympify_expression(li1.measure)*sympify_expression(li1.Gamma_factor)
    result_expoU_1 = sympify_expression(li1.exponent_U)
    result_expoF_1 = sympify_expression(li1.exponent_F)

    result_U_2 = sympify_expression(li2.U)
    result_F_2 = sympify_expression(li2.F)
    result_Nu_2 = sympify_expression(li2.numerator).subs('U',result_U_2).subs('F',result_F_2)\
                  *sympify_expression(li2.measure)*sympify_expression(li2.Gamma_factor)
    result_expoU_2 = sympify_expression(li2.exponent_U)
    result_expoF_2 = sympify_expression(li2.exponent_F)

    testcase.assertEqual( (result_U_1  - result_U_2 ).simplify() , 0 )
    testcase.assertEqual( (result_F_1  - result_F_2 ).simplify() , 0 )
    testcase.assertEqual( (result_Nu_1 - result_Nu_2).simplify() , 0 )
    testcase.assertEqual( (result_expoU_1 - result_expoU_2).simplify() , 0 )
    testcase.assertEqual( (result_expoF_1 - result_expoF_2).simplify() , 0 )


#@attr('active')
class TestPowerlist(unittest.TestCase):

    def tri1L_powers(self, power):
        loop_momenta = ['l']
        propagators = ['l**2', '(l-p1)**2', '(l+p2)**2']

        Feynman_parameters=['z1','z2','z3']

        rules = [('p1*p1','ssp1'),
                 ('p2*p2','ssp2'),
                 ('p1*p2','ssp3')]

        # compare against F, U, and Nu from SecDec3
        # For powers that are not directly supported by SecDec3, a second identical propagator was introduced
        # and the power split up, e.g. {-1/2+eps} -> {1/2+eps, -1}

        powerlist = sympify_expression([1,1,power])

        if powerlist[2].is_integer and powerlist[2].is_nonpositive:
            target_U = '''z1 + z2'''
            target_F = '''-(ssp1*z1*z2)'''
            target_integration_variables = sympify_expression(['z1', 'z2'])
        else:
            target_U = '''z1 + z2 + z3'''
            target_F = '''-((ssp1 + ssp2 + 2*ssp3)*z2*z3) - z1*(ssp1*z2 + ssp2*z3)'''
            target_integration_variables = sympify_expression(['z1', 'z2', 'z3'])

        symbols_U_F = sympify_expression(['U', 'F'])

        target_Nu = {  '1': '1'
                     , '0': '1'
                     , '-1':
                     '''3*ssp1*z1*z2 - 2*eps*ssp1*z1*z2 +
                     (z1 + z2)*(-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2)) -
                     eps*(z1 + z2)*(-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))'''
                     , '-2':
                     '''(-2 + eps)*(-10*ssp1**2*z1**2*z2**2 + 4*eps*ssp1**2*z1**2*z2**2 -
                     8*ssp1*z1*z2*(z1 + z2)*(-(ssp1*z2) - 2*ssp3*z2 -
                     ssp2*(z1 + z2)) + 4*eps*ssp1*z1*z2*(z1 + z2)*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2)) -
                     (z1 + z2)**2*(-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**2 +
                     eps*(z1 + z2)**2*(-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**2)'''
                     , '-3':
                     '''210*ssp1**3*z1**3*z2**3 - 214*eps*ssp1**3*z1**3*z2**3 +
                     72*eps**2*ssp1**3*z1**3*z2**3 - 8*eps**3*ssp1**3*z1**3*z2**3 +
                     270*ssp1**2*z1**2*z2**2*(z1 + z2)*(-(ssp1*z2) - 2*ssp3*z2 -
                     ssp2*(z1 + z2)) - 288*eps*ssp1**2*z1**2*z2**2*(z1 + z2)*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2)) +
                     102*eps**2*ssp1**2*z1**2*z2**2*(z1 + z2)*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2)) -
                     12*eps**3*ssp1**2*z1**2*z2**2*(z1 + z2)*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2)) +
                     90*ssp1*z1*z2*(z1 + z2)**2*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**2 -
                     111*eps*ssp1*z1*z2*(z1 + z2)**2*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**2 +
                     45*eps**2*ssp1*z1*z2*(z1 + z2)**2*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**2 -
                     6*eps**3*ssp1*z1*z2*(z1 + z2)**2*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**2 +
                     6*(z1 + z2)**3*(-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**
                     3 - 11*eps*(z1 + z2)**3*(-(ssp1*z2) - 2*ssp3*z2 -
                     ssp2*(z1 + z2))**3 + 6*eps**2*(z1 + z2)**3*
                     (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2))**3 -
                     eps**3*(z1 + z2)**3*(-(ssp1*z2) - 2*ssp3*z2 -
                     ssp2*(z1 + z2))**3'''
                     , '1+eps':
                     'z3**eps'
                     , 'eps':
                     'z3**(-1 + eps)'
                     , '-1+eps':
                     '''z3**(-1 + eps)*((-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*
                     (z1 + z2 + z3) - 2*eps*(-(ssp2*z1) -
                     (ssp1 + ssp2 + 2*ssp3)*z2)*(z1 + z2 + z3) -
                     3*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3)) +
                     3*eps*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3)))'''
                     , '-2+eps':
                     '''z3**(-1 + eps)*(2*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)**2*
                     (z1 + z2 + z3)**2 -
                     6*eps*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)**2*
                     (z1 + z2 + z3)**2 +
                     4*eps**2*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)**2*
                     (z1 + z2 + z3)**2 -
                     16*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*
                     (z1 + z2 + z3)*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3)) +
                     28*eps*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*
                     (z1 + z2 + z3)*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3)) -
                     12*eps**2*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*
                     (z1 + z2 + z3)*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3)) +
                     20*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3))**2 -
                     27*eps*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3))**2 +
                     9*eps**2*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3))**2)'''
                     , '1/2+eps':
                     '''z3**(1/2 + eps)
                     *1/2*((-1 - 4*eps)*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*(z1 + z2 + z3) +
                     3*(-1 + 2*eps)*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) - z1*(ssp1*z2 + ssp2*z3)))'''

                     , '-1/2+eps':
                     '''z3**(1/2 + eps)
                     *(((-1 + 4*eps)*(1 + 4*eps)*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)^2*
                     (z1 + z2 + z3)^2 - 2*(-1 + 4*eps)*(-5 + 6*eps)*
                     (-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*(z1 + z2 + z3)*
                     (-((ssp1 + ssp2 + 2*ssp3)*z2*z3) - z1*(ssp1*z2 + ssp2*z3)) +
                     (35 - 72*eps + 36*eps^2)*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) -
                     z1*(ssp1*z2 + ssp2*z3))^2)/4)
                     '''
                     ,'-3/2':
                       '''z3**(1/2)
                       *((-((-3 + 2*eps)*(-1 + 4*eps^2)*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)^3*
                       (z1 + z2 + z3)^3) + 3*(-7 + 4*eps)*(3 - 8*eps + 4*eps^2)*
                       (-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)^2*(z1 + z2 + z3)^2*
                       (-((ssp1 + ssp2 + 2*ssp3)*z2*z3) - z1*(ssp1*z2 + ssp2*z3)) -
                       3*(-3 + 2*eps)*(63 - 64*eps + 16*eps^2)*(-(ssp2*z1) - (ssp1 + ssp2 + 2*ssp3)*z2)*
                       (z1 + z2 + z3)*(-((ssp1 + ssp2 + 2*ssp3)*z2*z3) - z1*(ssp1*z2 + ssp2*z3))^2 +
                       (-693 + 956*eps - 432*eps^2 + 64*eps^3)*
                       (-((ssp1 + ssp2 + 2*ssp3)*z2*z3) - z1*(ssp1*z2 + ssp2*z3))^3)/8)
                       '''
                     }

        li = LoopIntegralFromPropagators(propagators, loop_momenta, powerlist=powerlist,
                                         replacement_rules=rules, Feynman_parameters=Feynman_parameters)

        number_of_derivatives = { '1': 0,
                                  '0': 0,
                                  '-1': 1,
                                  '-2': 2,
                                  '-3': 3,
                                  '1+eps': 0,
                                  'eps': 0,
                                  '-1+eps': 1,
                                  '-2+eps': 2,
                                  '1/2+eps': 1,
                                  '-1/2+eps': 2,
                                  '-3/2': 3 }

        # The powers *cannot* be compared against SecDec3 because the implementation is different!
        target_exponent_U = sum(powerlist) - number_of_derivatives[power] - (li.L + 1)*li.dimensionality/2
        target_exponent_F = - (sum(powerlist) + number_of_derivatives[power] - li.L*li.dimensionality/2)
        target_gamma = sp.gamma((sum(powerlist) - li.L*li.dimensionality/2)) * sp.exp(-sp.I*sp.pi * (sum(powerlist) - number_of_derivatives[power]))
        if (powerlist[2] + number_of_derivatives[power]) !=0:
            target_gamma /= sp.gamma(powerlist[2] + number_of_derivatives[power])

        result_U = sympify_expression(li.U)
        result_F = sympify_expression(li.F)
        result_Nu = sympify_expression(li.numerator).subs('U',result_U).subs('F',result_F)*sympify_expression(li.measure)

        # print "number_of_derivatives: ", li.number_of_derivatives
        # print "result_Nu = ", result_Nu
        # print "target_Nu = ", target_Nu[power]
        # print "ratio = ", (result_Nu/sympify_expression(target_Nu[power])).simplify()

        self.assertEqual( (result_U  - sympify_expression(target_U) ).simplify() , 0 )
        self.assertEqual( (result_F  - sympify_expression(target_F) ).simplify() , 0 )
        self.assertEqual( (result_Nu - sympify_expression(target_Nu[power])).simplify() , 0 )

        self.assertEqual( (li.exponent_U - target_exponent_U).simplify(), 0 )
        self.assertEqual( (li.exponent_F - target_exponent_F).simplify(), 0 )
        self.assertEqual( (li.Gamma_factor - target_gamma).simplify(), 0 )

        self.assertEqual( li.U.polysymbols, target_integration_variables)
        self.assertEqual( li.F.polysymbols, target_integration_variables)
        self.assertEqual( li.numerator.polysymbols, target_integration_variables + symbols_U_F)
        self.assertEqual( li.integration_variables, target_integration_variables)

    def test_one_power(self):
        self.tri1L_powers('1')

    def test_zero_power(self):
        self.tri1L_powers('0')

    #@attr('active')
    def test_negative_powers1(self):
        self.tri1L_powers('-1')

    def test_negative_powers2(self):
        self.tri1L_powers('-2')

    #@attr('active')
    @attr('slow')
    def test_negative_powers3(self):
        self.tri1L_powers('-3/2')

    @attr('slow')
    def test_negative_powers4(self):
        self.tri1L_powers('-3')

    def test_eps_power1(self):
        self.tri1L_powers('1+eps')

    def test_eps_power2(self):
        self.tri1L_powers('eps')

    def test_eps_power3(self):
        self.tri1L_powers('-1+eps')

    def test_eps_power4(self):
        self.tri1L_powers('-2+eps')

    #@attr('active')
    def test_eps_power5(self):
        self.tri1L_powers('1/2+eps')

    #@attr('active')
    def test_eps_power6(self):
        self.tri1L_powers('-1/2+eps')

    @attr('slow')
    def test_box_withnumerator2L(self):
        # SecDec3 -> loop/demos/4_box_withnumerator_2L
        loop_momenta = ['k1','k2']
        propagators = ['k1**2', '(k1+p2)**2', '(k1-p1)**2', '(k1-k2)**2', '(k2+p2)**2', '(k2-p1)**2',
                       '(k2+p2+p3)**2', '(k1+p3)**2']

        rules = [('p1*p1','0'),
                 ('p2*p2','0'),
                 ('p3*p3','0'),
                 ('p4*p4','0'),
                 ('p1*p2','s/2'),
                 ('p2*p3','t/2'),
                 ('p1*p3','-s/2-t/2')]

        powerlist = [1,1,1,1,1,1,1,-1]

        Feynman_parameters = ['z'+str(i) for i in range(1,9)]

        target_U = '''z4*(z5 + z6 + z7) + z1*(z4 + z5 + z6 + z7) +
        z2*(z4 + z5 + z6 + z7) + z3*(z4 + z5 + z6 + z7)'''

        target_F = '''-(s*z4*z5*z6) + s*z3*(-(z4*z5) - z5*z6) -
        z1*(s*z5*z6 + t*z4*z7) -
        s*z2*((z4 + z5)*z6 + z3*(z4 + z5 + z6 + z7))'''

        target_Nu = '''(-2*(z4*(z5 + z6 + z7) + z1*(z4 + z5 + z6 + z7) +
        z2*(z4 + z5 + z6 + z7) + z3*(z4 + z5 + z6 + z7))*
        (-(s*z5*z6) + z4*(s*z6 + t*(z5 + z6)) +
        t*z2*(z4 + z5 + z6 + z7) +
        z3*(s*(z4 + z5 + z6 + z7) +
        t*(z4 + z5 + z6 + z7))) -
        2*eps*(z4*(z5 + z6 + z7) + z1*(z4 + z5 + z6 + z7) +
        z2*(z4 + z5 + z6 + z7) + z3*(z4 + z5 + z6 + z7))*
        (-(s*z5*z6) + z4*(s*z6 + t*(z5 + z6)) +
        t*z2*(z4 + z5 + z6 + z7) +
        z3*(s*(z4 + z5 + z6 + z7) +
        t*(z4 + z5 + z6 + z7))) + 3*eps*(z4 + z5 + z6 + z7)*
        (-(s*z4*z5*z6) + s*z3*(-(z4*z5) - z5*z6) -
        z1*(s*z5*z6 + t*z4*z7) -
        s*z2*((z4 + z5)*z6 + z3*(z4 + z5 + z6 + z7))))'''

        target_gamma_factor = '-gamma(2+2*eps)'


        li = LoopIntegralFromPropagators(propagators, loop_momenta, powerlist=powerlist,
                                         replacement_rules=rules, Feynman_parameters=Feynman_parameters)

        result_U = sympify_expression(li.U)
        result_F = sympify_expression(li.F)
        result_Nu = sympify_expression(li.numerator).subs('U',result_U).subs('F',result_F)*sympify_expression(li.measure)
        result_gamma_factor = sympify_expression(li.Gamma_factor)

        self.assertEqual( (result_U  - sympify_expression(target_U) ).simplify() , 0 )
        self.assertEqual( (result_F  - sympify_expression(target_F) ).simplify() , 0 )
        self.assertEqual( (result_Nu - sympify_expression(target_Nu)).simplify() , 0 )
        self.assertEqual( (result_gamma_factor - sympify_expression(target_gamma_factor)).simplify() , 0 )

    @attr('slow')
    #@attr('active')
    def test_powerlist_with_mass(self):
        loop_momenta = ['k']
        propagators = ['k**2', '(k-p1)**2', '(k+p2)**2 - m**2']

        powerlist = [1,1,'-2+eps']

        rules = [('p1*p1','ssp1'),
                 ('p2*p2','ssp2'),
                 ('p1*p2','ssp3')]

        Feynman_parameters=['z1','z2','z3']

        li = LoopIntegralFromPropagators(propagators, loop_momenta, powerlist=powerlist,
                                         replacement_rules=rules, Feynman_parameters=Feynman_parameters)

        target_U = '''z1 + z2 + z3'''

        target_F = '''-(ssp1*z2*(z1 + z3)) + z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3))'''

        target_Nu = '''z3**(-1 + eps)*(2*(z1 + z2 + z3)**2*
        (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2) + m**2*z3 +
        m**2*(z1 + z2 + z3))**2 - 6*eps*(z1 + z2 + z3)**2*
        (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2) + m**2*z3 +
        m**2*(z1 + z2 + z3))**2 + 4*eps**2*(z1 + z2 + z3)**2*
        (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2) + m**2*z3 +
        m**2*(z1 + z2 + z3))**2 + 4*m**2*(z1 + z2 + z3)**2*
        (-(ssp1*z2*(z1 + z3)) + z3*(-2*ssp3*z2 -
        ssp2*(z1 + z2) + m**2*(z1 + z2 + z3))) -
        4*eps*m**2*(z1 + z2 + z3)**2*(-(ssp1*z2*(z1 + z3)) +
        z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3))) - 16*(z1 + z2 + z3)*
        (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2) + m**2*z3 +
        m**2*(z1 + z2 + z3))*(-(ssp1*z2*(z1 + z3)) +
        z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3))) + 28*eps*(z1 + z2 + z3)*
        (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2) + m**2*z3 +
        m**2*(z1 + z2 + z3))*(-(ssp1*z2*(z1 + z3)) +
        z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3))) - 12*eps**2*(z1 + z2 + z3)*
        (-(ssp1*z2) - 2*ssp3*z2 - ssp2*(z1 + z2) + m**2*z3 +
        m**2*(z1 + z2 + z3))*(-(ssp1*z2*(z1 + z3)) +
        z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3))) +
        20*(-(ssp1*z2*(z1 + z3)) + z3*(-2*ssp3*z2 -
        ssp2*(z1 + z2) + m**2*(z1 + z2 + z3)))**2 -
        27*eps*(-(ssp1*z2*(z1 + z3)) +
        z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3)))**2 +
        9*eps**2*(-(ssp1*z2*(z1 + z3)) +
        z3*(-2*ssp3*z2 - ssp2*(z1 + z2) +
        m**2*(z1 + z2 + z3)))**2)'''

        target_gamma_factor = 'gamma(-2+2*eps)/gamma(eps) * exp(-I*pi * (eps+2))'

        result_U = sympify_expression(li.U)
        result_F = sympify_expression(li.F)
        result_Nu = sympify_expression(li.numerator).subs('U',result_U).subs('F',result_F)*sympify_expression(li.measure)
        result_gamma_factor = sympify_expression(li.Gamma_factor)

        self.assertEqual( (result_U  - sympify_expression(target_U) ).simplify() , 0 )
        self.assertEqual( (result_F  - sympify_expression(target_F) ).simplify() , 0 )
        self.assertEqual( (result_Nu - sympify_expression(target_Nu)).simplify() , 0 )
        self.assertEqual( (result_gamma_factor - sympify_expression(target_gamma_factor)).simplify() , 0 )


    def test_powerlist_cutconstruct(self):
        # like SecDec3 -> loop/demos/2_triangle_2L, but with non-trivial powerlist
        internal_lines = [['m',[3,4]],['m',[4,5]],['m',[3,5]],[0,[1,2]],[0,[4,1]],[0,[2,5]]]
        external_lines = [['p1',1], ['p2',2], ['-p1',3], ['-p2',3]]
        powerlist = [1,2,3,1,4,1]

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p1*p2','s/2'),
                  ('m**2','ms1')]

        Feynman_parameters=['z' + str(i) for i in range(1,7)]

        li = LoopIntegralFromGraph(internal_lines, external_lines, powerlist=powerlist,
                                   replacement_rules=rules, Feynman_parameters=Feynman_parameters)

        target_U = '''z1*z2 + z2*z3 + z1*z4 + z2*z4 + z3*z4 + z1*z5 +
        z2*z5 + z3*z5 + z1*z6 + z2*z6 + z3*z6'''

        target_F = '''-(s*z1*z2*z3) - s*z1*z3*z4 - s*z1*z3*z5 -
        s*z2*z3*z5 - s*z1*z2*z6 - s*z1*z3*z6 -
        s*z1*z5*z6 - s*z2*z5*z6 - s*z3*z5*z6 +
        ms1*(z1 + z2 + z3)*(z3*(z4 + z5 + z6) +
        z1*(z2 + z4 + z5 + z6) + z2*(z3 + z4 + z5 + z6))'''

        target_Nu = '''z2*z3**2*z5**3'''

        target_gamma = sp.gamma('8 + 2*eps')/12

        result_U = sympify_expression(li.U)
        result_F = sympify_expression(li.F)
        result_Nu = sympify_expression(li.numerator)*sympify_expression(li.measure)
        result_gamma = li.Gamma_factor

        self.assertEqual( li.external_momenta, sympify_expression(['p1','p2']) )
        self.assertEqual( (result_U  - sympify_expression(target_U) ).simplify() , 0 )
        self.assertEqual( (result_F  - sympify_expression(target_F) ).simplify() , 0 )
        self.assertEqual( (result_Nu - sympify_expression(target_Nu)).simplify() , 0 )
        self.assertEqual( (result_gamma  - target_gamma ).simplify() , 0 )

    #@attr('active')
    def test_multiple_negative_powers(self):
        internal_lines = [[0,[1,2]], [0,[2,3]], [0,[3,4]], [0,[4,1]]]
        external_lines = [['p1',1], ['p2',2], ['p3',3], ['-p1-p2-p3',4]]
        powerlist = [-1,1,-1,1]

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p3*p3','0'),
                  ('p1*p2','s/2'),
                  ('p2*p3','t/2'),
                  ('p1*p3','-s/2-t/2')]

        Feynman_parameters=['dummy1','z1','dummy2','z2']

        li = LoopIntegralFromGraph(internal_lines, external_lines, powerlist=powerlist,
                                   replacement_rules=rules, Feynman_parameters=Feynman_parameters)

        target_U = '''z1 + z2'''

        target_F = '''-(s*z1*z2)'''

        target_Nu = '''20*s**2*z1**2*z2**2 - 18*eps*s**2*z1**2*z2**2 +
        4*eps**2*s**2*z1**2*z2**2 + 2*s*t*z1*z2*(z1 + z2)**2 -
        eps*s*t*z1*z2*(z1 + z2)**2'''

        target_gamma = sp.gamma('eps - 2')

        target_measure = 1

        target_measure_symbols = sympify_expression([ 'z1', 'z2', 'U', 'F' ])

        result_U = sympify_expression(li.U)
        result_F = sympify_expression(li.F)
        result_Nu = sympify_expression(li.numerator).subs('U',result_U).subs('F',result_F)*sympify_expression(li.measure)
        result_gamma = li.Gamma_factor
        result_measure = sympify_expression(li.measure)
        result_measure_symbols = li.measure.symbols

        self.assertEqual( (result_U  - sympify_expression(target_U) ).simplify() , 0 )
        self.assertEqual( (result_F  - sympify_expression(target_F) ).simplify() , 0 )
        self.assertEqual( (result_Nu - sympify_expression(target_Nu)).simplify() , 0 )
        self.assertEqual( (result_gamma  - target_gamma ).simplify() , 0 )
        self.assertEqual( (result_measure - sympify_expression(target_measure)).simplify() , 0 )
        self.assertEqual( result_measure_symbols, target_measure_symbols)

    #@attr('active')
    def test_one_zero_and_one_negative_power(self):
        external_lines = [['p1',1], ['p2',2], ['p3',3], ['-p1-p2-p3',4]]

        internal_lines_original = [[0,[1,2]], [0,[2,3]], [0,[3,4]], [0,[4,1]]]
        powerlist_original = [-1,3,0,1]

        internal_lines_swapped = [[0,[3,4]], [0,[2,3]], [0,[1,2]], [0,[4,1]]]
        powerlist_swapped = [0,3,-1,1]

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p3*p3','0'),
                  ('p1*p2','s/2'),
                  ('p2*p3','t/2'),
                  ('p1*p3','-s/2-t/2')]

        Feynman_parameters=['dummy1','z1','dummy2','z2']

        for internal_lines,powerlist in zip([internal_lines_original,internal_lines_swapped],
                                            [     powerlist_original,     powerlist_swapped]):

            li = LoopIntegralFromGraph(internal_lines, external_lines, powerlist=powerlist,
                                       replacement_rules=rules, Feynman_parameters=Feynman_parameters)

            target_U = '''z1 + z2'''

            target_F = '''-(s*z1*z2)'''

            target_Nu = '''z1**2*(s*z1*z2 - 2*eps*s*z1*z2)'''

            target_measure = '''z1**2'''

            target_gamma = sp.gamma('1 + eps')/2

            result_U = sympify_expression(li.U)
            result_F = sympify_expression(li.F)
            result_Nu = sympify_expression(li.numerator).subs('U',result_U).subs('F',result_F)*sympify_expression(li.measure)
            result_gamma = li.Gamma_factor
            result_measure = sympify_expression(li.measure)
            result_measure_symbols = li.measure.symbols
            target_measure_symbols = sympify_expression([ 'z1', 'z2', 'U', 'F' ])

            self.assertEqual( (result_U  - sympify_expression(target_U) ).simplify() , 0 )
            self.assertEqual( (result_F  - sympify_expression(target_F) ).simplify() , 0 )
            self.assertEqual( (result_Nu - sympify_expression(target_Nu)).simplify() , 0 )
            self.assertEqual( (result_gamma  - target_gamma ).simplify() , 0 )
            self.assertEqual( (result_measure - sympify_expression(target_measure)).simplify() , 0 )
            self.assertEqual( result_measure_symbols, target_measure_symbols)

    @attr('slow')
    #@attr('active')
    def test_compare_tensor_integral_to_inverse_propagator(self):
        loop_momenta = ['l']
        propagators = ['l**2', '(l-p1)**2', '(l+p2)**2', '(l+p2+p3)**2','l*l']

        powerlist1 = [0,1,1,0,0]
        numerator1 = 'l(mu)*l(mu) * l(nu)*l(nu) * l(rho)*l(rho)'
        indices1 = ['mu','nu','rho']

        powerlist2 = [-3,1,1,0,0]
        numerator2 = '1'
        indices2 = []

        powerlist3 = [0,1,1,0,-3]
        numerator3 = '1'
        indices3 = []

        external_momenta = ['p1', 'p2']

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p3*p3','0'),
                  ('p1*p2','ssp1/2'),
                  ('p2*p3','ssp2/2'),
                  ('p1*p3','-ssp1/2-ssp2/2')]

        li1 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        li3 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator3,
                                          Lorentz_indices=indices3, powerlist=powerlist3, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)
        compare_two_loop_integrals(self, li1, li3)
        compare_two_loop_integrals(self, li2, li3)

    @attr('slow')
    #@attr('active')
    def test_compare_tensor_integral_to_inverse_propagator2(self):
        loop_momenta = ['l']
        propagators = ['l**2', '(l-p1)**2', '(l+p2)**2']

        powerlist1 = [0,1,1]
        numerator1 = 'l(mu)*l(mu) * l(nu)*l(nu) * l(rho)*l(rho) * l(sigma)*l(sigma) '
        indices1 = ['mu','nu','rho','sigma']

        powerlist2 = [-4,1,1]
        numerator2 = '1'
        indices2 = []

        external_momenta = ['p1', 'p2']

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p3*p3','0'),
                  ('p1*p2','ssp1/2'),
                  ('p2*p3','ssp2/2'),
                  ('p1*p3','-ssp1/2-ssp2/2')]

        li1 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)

    @attr('slow')
    #@attr('active')
    def test_combined_tensor_integral_and_inverse_propagator_1(self):
        loop_momenta = ['l']
        propagators = ['l**2', '(l-p1)**2', '(l+p2)**2', '(l+p2+p3)**2']

        powerlist1 = [1,1,1,0]
        numerator1 = '(l(a)+p2(a)+p3(a))*(l(a)+p2(a)+p3(a)) * l(mu)*p2(mu) * l(nu)'
        indices1 = ['mu', 'nu', 'a']

        powerlist2 = [1,1,1,-1]
        numerator2 = 'l(mu)*p2(mu) * l(nu)'
        indices2 = ['mu', 'nu']

        external_momenta = ['p1', 'p2', 'p3']

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p3*p3','0'),
                  ('p1*p2','ssp1/2'),
                  ('p2*p3','ssp2/2'),
                  ('p1*p3','-ssp1/2-ssp2/2')]

        li1 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)


    #@attr('active')
    def test_combined_tensor_integral_and_inverse_propagator_2(self):
        loop_momenta = ['l']
        propagators = ['l**2', '(l-p1)**2', '(l+p2)**2']

        powerlist1 = [0,1,1]
        numerator1 = 'l(mu)*l(mu) * l(nu)*l(nu) * l(rho)'
        indices1 = ['mu','nu','rho']

        powerlist2 = [-2,1,1]
        numerator2 = 'l(rho)'
        indices2 = ['rho']

        external_momenta = ['p1', 'p2']

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p1*p2','ssp1/2')]

        li1 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)


    #@attr('active')
    def test_combined_tensor_integral_and_inverse_propagator_3(self):
        loop_momenta = ['l']
        propagators = ['l**2', '(l-p1)**2', '(l+p2)**2']

        powerlist1 = [1,'1/2+eps',1]
        numerator1 = '(l(mu)-p1(mu))*(l(mu)-p1(mu)) * l(nu)'
        indices1 = ['mu','nu']

        powerlist2 = [1,'-1/2+eps',1]
        numerator2 = 'l(nu)'
        indices2 = ['nu']

        external_momenta = ['p1', 'p2']

        rules = [ ('p1*p1','0'),
                  ('p2*p2','0'),
                  ('p1*p2','ssp1/2')]

        li1 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)


    @attr('slow')
    #@attr('active')
    def test_combined_tensor_integral_and_inverse_propagator_4(self):
        loop_momenta = ['l', 'k']

        propagators1 = ['l**2', '(l+p1)**2', 'k**2', '(k-p1)**2', '(l+k)**2']
        powerlist1 = [1,-1,1,1,1]
        numerator1 = '(l(mu) - k(mu)) * (l(mu) - k(mu))'
        indices1 = ['mu']

        propagators2 = ['l**2', '(l+p1)**2', 'k**2', '(k-p1)**2', '(l+k)**2', '(l-k)**2']
        powerlist2 = [1,-1,1,1,1,-1]
        numerator2 = '1'
        indices2 = []

        external_momenta = ['p1']

        rules = [ ('p1*p1','p1sqr')]

        li1 = LoopIntegralFromPropagators(propagators1, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators2, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)


    #@attr('active')
    def test_inverse_linear_propagator(self):
        loop_momenta = ['l']

        propagators1 = ['l**2', '(l+p1)**2']
        powerlist1 = [1,1]
        numerator1 = 'l(mu)*p1(mu) * l(nu)*p1(nu) * l(rho)*p1(rho)'
        indices1 = ['mu', 'nu', 'rho']

        propagators2 = ['l**2', '(l+p1)**2', 'l*p1']
        powerlist2 = [1,1,-3]
        numerator2 = '1'
        indices2 = []

        external_momenta = ['p1']

        rules = [ ('p1*p1','p1sqr')]

        li1 = LoopIntegralFromPropagators(propagators1, loop_momenta, external_momenta, numerator=numerator1,
                                          Lorentz_indices=indices1, powerlist=powerlist1,
                                          replacement_rules=rules)

        li2 = LoopIntegralFromPropagators(propagators2, loop_momenta, external_momenta, numerator=numerator2,
                                          Lorentz_indices=indices2, powerlist=powerlist2, replacement_rules=rules)

        compare_two_loop_integrals(self, li1, li2)


    #@attr('active')
    # this integral is known analytically -> candidate to check numerical result
    def test_Y_integral(self):

        loop_momenta = ['k1', 'k2']
        propagators = ['k1**2', 'k2**2', '(k1+k2)**2', '-b2*k2-1', '-b1*k1-1', 'b2*k1-2',
                       '(-aij/2-1/(2*aij))*b2*(2*k2+k1)-b1*(2*k2+k1)']
        powerlist = [1,1,1,1,1,1,-1]

        rules = [('b1*b1','1'),
                 ('b2*b2','1'),
                 ('b1*b2','-aij/2-1/(2*aij)')]

        Feynman_parameters=['z' + str(i) for i in range(1,8)]

        li = LoopIntegralFromPropagators(propagators, loop_momenta, powerlist=powerlist,
                                         replacement_rules=rules, Feynman_parameters=Feynman_parameters)


        target_U = sympify_expression('''z2*z3 + z1*(z2 + z3)''')

        target_F = sympify_expression('''(4*aij*z5*(z2*z6 + z3*(z4 + z6)) +
                    4*aij**3*z5*(z2*z6 + z3*(z4 + z6)) +
                    2*aij**2*(z3*(2*z4**2 + 2*z5**2 + 4*z4*z6 + 2*z6**2) +
                    2*z1*(z4**2 + 4*z2*(z4 + z5 + 2*z6) +
                    4*z3*(z4 + z5 + 2*z6)) +
                    z2*(2*z5**2 + 2*z6**2 + 8*z3*(z4 + z5 + 2*z6))))/(16*aij**2)''')

        target_Nu = sympify_expression('''-((z2*z3 + z1*(z2 + z3))*(-2*z2*z5 + 2*z3*z5 +
                    2*aij**2*(2*z2*z5 - 2*z3*z5) -
                    aij**4*(2*z2*z5 - 2*z3*z5)))/(16*aij**2) -
                    (eps*(z2*z3 + z1*(z2 + z3))*(-2*z2*z5 + 2*z3*z5 +
                    2*aij**2*(2*z2*z5 - 2*z3*z5) -
                    aij**4*(2*z2*z5 - 2*z3*z5)))/(8*aij**2)''')

        result_U = sympify_expression(li.U)
        result_F = sympify_expression(li.F)
        result_Nu = sympify_expression(li.numerator).subs('U',result_U).subs('F',result_F)*sympify_expression(li.measure)

        # The `target_Nu` produced with SecDec3 has an additional factor of U in the numerator...
        result_Nu *= result_U

        self.assertEqual( (result_U  - target_U ).simplify() , 0 )
        self.assertEqual( (result_F  - target_F ).simplify() , 0 )
        self.assertEqual( (result_Nu - target_Nu).simplify() , 0 )

    #@attr('active')
    def test_powerlist_regulators_bubble_1L_from_graph(self):
        internal_lines = [['mt',[1,2]],['mt',[2,1]]]
        external_lines = [['p1',1],['p2',2]]
        powerlist=['1+n1','1+n2']
        regulators=['eps','n1','n2']
        replacement_rules = [
                                ('p1*p1', 'msq'),
                                ('p2*p2', 'msq'),
                                ('mt**2', 'mtsq'),
                            ]
        result_U = "x0+x1"
        result_exponent_U = "2*eps+n1+n2-2"
        result_F = "mtsq*(x0**2+2*x0*x1+x1**2)-x0*x1*msq"
        result_exponent_F = "-eps-n1-n2"
        result_measure = "x0**n1*x1**n2"
        loop_integral = LoopIntegralFromGraph(internal_lines=internal_lines, external_lines=external_lines, powerlist=powerlist, regulators=regulators, replacement_rules=replacement_rules)
        # uf_from_graph_generic(self, internal_lines, external_lines, 1, result_U, result_F, regulators=regulators, powerlist=powerlist)

        self.assertEqual(loop_integral.L,1)
        self.assertEqual((sympify_expression(loop_integral.U)-sympify_expression(result_U)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.exponent_U)-sympify_expression(result_exponent_U)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.F)-sympify_expression(result_F)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.exponent_F)-sympify_expression(result_exponent_F)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.measure)-sympify_expression(result_measure)).simplify(),0)

    #@attr('active')
    def test_powerlist_regulators_bubble_1L_from_propagators(self):
        propagators = ["(k**2-mtsq)", "((p1-k)**2-mtsq)"]
        powerlist=['1+n1','1+n2']
        regulators=['eps','n1','n2']
        replacement_rules = [
                                ('p1*p1', 'msq'),
                                ('p2*p2', 'msq'),
                                ('mt**2', 'mtsq'),
                            ]
        result_U = "x0+x1"
        result_exponent_U = "2*eps+n1+n2-2"
        result_F = "mtsq*(x0**2+2*x0*x1+x1**2)-x0*x1*msq"
        result_exponent_F = "-eps-n1-n2"
        result_measure = "x0**n1*x1**n2"
        loop_integral = LoopIntegralFromPropagators(propagators=propagators, loop_momenta=["k"], powerlist=powerlist, regulators=regulators, replacement_rules=replacement_rules)
        # uf_from_graph_generic(self, internal_lines, external_lines, 1, result_U, result_F, regulators=regulators, powerlist=powerlist)

        self.assertEqual(loop_integral.L,1)
        self.assertEqual((sympify_expression(loop_integral.U)-sympify_expression(result_U)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.exponent_U)-sympify_expression(result_exponent_U)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.F)-sympify_expression(result_F)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.exponent_F)-sympify_expression(result_exponent_F)).simplify(),0)
        self.assertEqual((sympify_expression(loop_integral.measure)-sympify_expression(result_measure)).simplify(),0)

#@attr('active')
class TestUF_LoopPackageFromPropagators(unittest.TestCase):
    def test_loop_package_twice(self):
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp(prefix="pysecdec")
        try:
            os.chdir(tmpdir)
            li = LoopIntegralFromPropagators(
                propagators=["k**2", "(q - k)**2"],
                loop_momenta=["k"],
                external_momenta=["q"],
                replacement_rules=[("q**2", "1")]
            )
            loop_package("bubble1", li, [0])
            loop_package("bubble2", li, [0])
            # The test is that no exceptions are thrown up to
            # this point.
        finally:
            os.chdir(cwd)
            shutil.rmtree(tmpdir)

class TestUF_LoopPackageFromGraph(unittest.TestCase):
    def test_loop_package_twice(self):
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp(prefix="pysecdec")
        try:
            os.chdir(tmpdir)
            li = LoopIntegralFromGraph(
                internal_lines=[[0, [1, 2]], [0, [1,2]]],
                external_lines=[['p', 1], ['p', 2]],
                replacement_rules=[('p*p', '-1')]
            )
            loop_package("bubble1", li, [0])
            loop_package("bubble2", li, [0])
            # The test is that no exceptions are thrown up to
            # this point.
        finally:
            os.chdir(cwd)
            shutil.rmtree(tmpdir)
