"""Unit tests for the U, F routines"""

from .loop_integral import *
import sympy as sp
import unittest
from nose.plugins.attrib import attr

def uf_from_propagators_generic(test_case, loop_momenta, propagators, result_u, result_f):
    loop_integral = LoopIntegral_from_propagators(propagators, loop_momenta, Feynman_parameters='Feynman')

    u,f = loop_integral.U, loop_integral.F

    sympy_u = sp.sympify(str(u))
    sympy_f = sp.sympify(str(f))

    result_u = sp.sympify(str(result_u).replace('x','Feynman'))
    result_f = sp.sympify(str(result_f).replace('x','Feynman'))

    zerou = (sympy_u - result_u).expand()
    zerof = (sympy_f - result_f).expand()

    test_case.assertEqual(zerou,0)
    test_case.assertEqual(zerof,0)

    expo_u, expo_f = loop_integral.exponent_U, loop_integral.exponent_F
    test_case.assertEqual(expo_u, len(propagators) - sp.sympify('2-eps') * (1 + len(loop_momenta)))
    test_case.assertEqual(expo_f, -(len(propagators) - sp.sympify('2-eps') * len(loop_momenta)))

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
                                LoopIntegral_from_propagators, ['k1**2'], ['k1'], Feynman_parameters=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*loop_momenta.*symbol', \
                                LoopIntegral_from_propagators, ['k1**2', '(k1+k2)**2'], ['k1', '2*k2'], Feynman_parameters=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*external_momenta.*symbol', \
                                LoopIntegral_from_propagators, ['k1**2'], ['k1', 'k2'], ['alpha*p1'], Feynman_parameters=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*Lorentz_indices.*symbol.*number', \
                                LoopIntegral_from_propagators, ['k1**2'], ['k1', 'k2'], Lorentz_indices=['alpha*delta'])
        self.assertRaisesRegexp(AssertionError, '.*metric_tensor.*symbol', \
                                LoopIntegral_from_propagators, ['k1**2'], ['k1', 'k2'], metric_tensor='4')
        self.assertRaisesRegexp(AssertionError, '.*propagators.*at most quadratic', \
                                LoopIntegral_from_propagators, ['k1**2', 'p1*(k1+k2)**2'], ['k1','k2'], ['p1','p2'])
        self.assertRaisesRegexp(AssertionError, '.*replacement_rules.*list of tuples', \
                                LoopIntegral_from_propagators, ['k1**2', '(k1+k2)**2'], ['k1','k2'], replacement_rules=['z0','z1'])
        self.assertRaisesRegexp(AssertionError, '.*replacement_rules.*list of tuples', \
                                LoopIntegral_from_propagators, ['k1**2', '(k1+k2)**2'], ['k1','k2'], replacement_rules=[('k1*k1', 0,'z1')])
        self.assertRaisesRegexp(AssertionError, '.*replacement_rules.*at most quadratic', \
                                LoopIntegral_from_propagators, ['k1**2', '(k1+k2)**2'], ['k1','k2'], replacement_rules=[('k1*k1*k2', 0)])

    #@attr('active')
    def test_replacement_rules_box_2L(self):
        k1, k2, p1, p2, p3, p4, s, t = sp.symbols('k1 k2 p1 p2 p3 p4 s t')

        loop_momenta = [k1, k2]
        propagators = [k1**2,(k1+p2)**2,(k1-p1)**2,(k1-k2)**2,(k2+p2)**2,(k2-p1)**2,(k2+p2+p3)**2]
        replacement_rules = [(p1*p1, 0),
                             (p2*p2, 0),
                             (p3*p3, 0),
                             (p4*p4, 0),
                             (p1*p2, s/2),
                             (p2*p3, t/2),
                             (p1*p3, -s/2-t/2)]
        loop_integral = LoopIntegral_from_propagators(propagators, loop_momenta, replacement_rules=replacement_rules, Feynman_parameters='z')
        U = sp.sympify(loop_integral.U)
        F = sp.sympify(loop_integral.F)

        target_U = sp.sympify('z3*(z4 + z5 + z6) + z0*(z3 + z4 + z5 + z6) + z1*(z3 + z4 + z5 + z6) + z2*(z3 + z4 + z5 + z6)')
        target_F = sp.sympify('-(s*z3*z4*z5) + s*z2*(-(z3*z4) - z4*z5) - z0*(s*z4*z5 + t*z3*z6) - s*z1*((z3 + z4)*z5 + z2*(z3 + z4 + z5 + z6))')

        self.assertEqual( (target_U - U).simplify() , 0)
        self.assertEqual( (target_F - F).simplify() , 0)

        L = 2
        D = sp.sympify('4-2*eps')
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
        loop_integral = LoopIntegral_from_propagators(props, loop_momenta, replacement_rules=replacement_rules, Feynman_parameters=['z1','z2','z3','z4'])
        U = sp.sympify(loop_integral.U)
        F = sp.sympify(loop_integral.F)

        target_U = sp.sympify('z1 + z2 + z3 + z4')
        target_F = sp.sympify('-(t*z1*z3) - (s1*z1 + s*z2)*z4 + m**2*z1*(z1 + z2 + z3 + z4)')

        self.assertEqual( (target_U - U).simplify() , 0)
        self.assertEqual( (target_F - F).simplify() , 0)

class TestNumerator(unittest.TestCase):
    #@attr('active')
    def test_double_index_notation(self):
        mu, nu = sp.symbols('mu nu')
        numerator = 'k1(mu)*(k2(mu) + p1(mu)) + k2(mu)*k2(nu)*p1(mu)*p2(nu)'
        indices = ['mu', 'nu']
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator, Lorentz_indices=indices)
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

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator, Lorentz_indices=[mu])
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

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator, Lorentz_indices=indices)
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [[(0,mu),(1,mu),(1,1),(1,2)]]

        self.assertEqual(len(tensors), 1)
        for i in range(1):
            print(i)
            self.assertEqual(tensors[i], target_tensors[i])

    #@attr('active')
    def test_tensor_numerator_bubble_1L(self):
        numerator = sp.sympify('k(1)*k(2)*k(3)*k(4)')
        contracted_numerator = numerator * sp.sympify('p(1)*p(2)*p(3)*p(4)*const')
        partially_contracted_numerator = numerator * sp.sympify('p(3)*p(4)')
        indices = [1, 2, 3, 4]
        loop_momenta = ['k']
        external_momenta = ['p']
        propagators = ['k**2', '(k - p)**2']

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                           numerator=numerator, Feynman_parameters=['x1','x2'],
                                           Lorentz_indices=indices, dimensionality='D')
        li_contracted = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                                      numerator=contracted_numerator, dimensionality='D',
                                                      Feynman_parameters=['x1','x2'], Lorentz_indices=indices)
        li_partially_contracted = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                                                numerator=partially_contracted_numerator,
                                                                Feynman_parameters=['x1','x2'], Lorentz_indices=indices,
                                                                replacement_rules=[('p*p', 'm**2')], dimensionality='D')
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [[(0,1),(0,2),(0,3),(0,4)]]

        self.assertEqual(tensors, target_tensors)

        numerator = sp.sympify(li.numerator) * li.Gamma_factor
        contracted_numerator = sp.sympify(li_contracted.numerator) * li.Gamma_factor
        partially_contracted_numerator = sp.sympify(li_partially_contracted.numerator) * li.Gamma_factor
        scalar_factor_insertation = { # N_nu = 2, L = 1 in this example
                                          'scalar_factor(0)': '1/(-2)**(0/2)*(2 - D*1/2 - 2/2)*(2 - D*1/2 - 4/2)*gamma(2 - D*1/2 - 4/2)*F**(0/2)',
                                          'scalar_factor(2)': '1/(-2)**(2/2)*(2 - D*1/2 - 4/2)*gamma(2 - D*1/2 - 4/2)*F**(2/2)',
                                          'scalar_factor(4)': '1/(-2)**(4/2)*gamma(2 - D*1/2 - 4/2)*F**(4/2)'
                                    }
        target_numerator = sp.sympify('''
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
        target_contracted_numerator = sp.sympify('''
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
                                                 ''').subs(scalar_factor_insertation) * sp.sympify('const')
        target_partially_contracted_numerator = sp.sympify('''
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
        numerator = sp.sympify('k(1)*k(2)*k(3)*k(4)*p(1)*p(2)*p(3)*p(4)')
        indices = [1, 2, 3, 4]
        loop_momenta = ['k']
        external_momenta = ['p']
        propagators = ['k**2', '(k - p)**2']
        replacement_rules = [('p**2', 'm**2')]

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator, Feynman_parameters=['x1','x2'], replacement_rules=replacement_rules, Lorentz_indices=indices)

        target_numerator = sp.sympify('''
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
                                                      'scalar_factor(0)': '1/(-2)**(0/2)*(2 - (4-2*eps)*1/2 - 2/2)*(2 - (4-2*eps)*1/2 - 4/2)*gamma(2 - (4-2*eps)*1/2 - 4/2)*F**(0/2)',
                                                      'scalar_factor(2)': '1/(-2)**(2/2)*(2 - (4-2*eps)*1/2 - 4/2)*gamma(2 - (4-2*eps)*1/2 - 4/2)*F**(2/2)',
                                                      'scalar_factor(4)': '1/(-2)**(4/2)*gamma(2 - (4-2*eps)*1/2 - 4/2)*F**(4/2)'
                                                })

        self.assertEqual( (sp.sympify(li.numerator)*li.Gamma_factor - target_numerator).simplify() , 0 )

    #@attr('active')
    def test_2L_box_with_numerator(self):
        propagators = ['k1**2','(k1+p2)**2','(k1-p1)**2','(k1-k2)**2','(k2+p2)**2','(k2-p1)**2','(k2+p2+p3)**2']
        numerator = '2*k1(mu)*p3(mu) + eps'
        loop_momenta = ['k1','k2']
        external_momenta = ['p1','p2','p3','p4']

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator, Feynman_parameters=['z%i'%i for i in range(1,7+1)], Lorentz_indices=['mu'])

        # comparison with Mathematica implementation of SecDec
        target_U = sp.sympify('''
                                   z4*(z5 + z6 + z7) + z1*(z4 + z5 + z6 + z7) +
                                   z2*(z4 + z5 + z6 + z7) + z3*(z4 + z5 + z6 + z7)
                              ''')
        target_F = sp.sympify('''
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
        target_numerator = sp.sympify('''
                                               2*(-(p3(mu)* p3(mu)*z4*z7)
                                             - p2(mu)* p3(mu)*(z4*(z5 + z7)
                                             + z2*(z4 + z5 + z6 + z7))
                                             + p1(mu)* p3(mu)*(z4*z6 + z3*(z4 + z5 + z6 + z7)))
                                             + eps * U
                                      ''') * sp.sympify('gamma(7 - (4-2*eps)*2/2)')
        # SecDec puts the ``(-1)**N_nu`` into the prefactor while `LoopIntegral` puts it into the numerator.
        target_numerator = - target_numerator

        self.assertEqual( (sp.sympify(li.U) - target_U).simplify() , 0 )
        self.assertEqual( (sp.sympify(li.F) - target_F  ).simplify() , 0 )
        self.assertEqual( (sp.sympify(li.numerator) * li.Gamma_factor - target_numerator).simplify() , 0 )

        replacement_rules = np.array([('p1**2', 0),
                                      ('p2**2', 0),
                                      ('p3**2', 0),
                                      ('p4**2', 0),
                                      ('p1*p2', 's/2'),
                                      ('p2*p3', 't/2'),
                                      ('p1*p3', '-s/2 - t/2')])

        li_with_replacement_rules = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                                                  numerator=numerator, Feynman_parameters=['z%i'%i for i in range(1,7+1)],
                                                                  replacement_rules=replacement_rules, Lorentz_indices=['mu'])

        target_numerator_with_replacements = sp.sympify('''
                                                             2*(
                                                                    - (t/2)*(z4*(z5 + z7)
                                                                    + z2*(z4 + z5 + z6 + z7))
                                                                    + (-s/2-t/2)*(z4*z6 + z3*(z4 + z5 + z6 + z7))
                                                               ) + eps * U
                                                        ''') * sp.sympify('gamma(7 - (4-2*eps)*2/2)')
        # SecDec puts the ``(-1)**N_nu`` into the prefactor while `LoopIntegral` puts it into the numerator.
        target_numerator_with_replacements = - target_numerator_with_replacements

        self.assertEqual( (sp.sympify(li_with_replacement_rules.numerator) * li.Gamma_factor - target_numerator_with_replacements).simplify() , 0 )

    #@attr('active')
    @attr('slow')
    def test_rank3_numerator_2L_double_triangle(self):
        k1, k2, p1, p2, p3, p4, s, t, mu, nu, scalar_factor, eps, F = sp.symbols('k1 k2 p1 p2 p3 p4 s t mu nu scalar_factor eps F')
        sp11, sp12, sp13, sp22, sp23, sp33 = sp.symbols('sp11 sp12 sp13 sp22 sp23 sp33')

        z = sp.sympify(['z%i'%i for i in range(6+1)]) # Feynman parameters (only use z[1], z[2], ..., z[6] but not z[0])
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
        numerator = k1(mu)*k2(nu)*k1(mu)*p2(nu)*2

        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                           numerator=numerator, Feynman_parameters=z[1:],
                                           Lorentz_indices=indices,
                                           replacement_rules=replacement_rules)

        # need to insert the `scalar_factor` (and the `F` therein) when comparing with SecDec
        # `scalar_factor` = `1/(-2)**(r/2)*Gamma(N_nu - dim*L/2 - r/2)*F**(r/2)
        N_nu = sp.sympify(len(propagators))
        dim = sp.sympify(4 - 2*eps)
        L = sp.sympify(2)
        numerator = sp.sympify(li.numerator) * li.Gamma_factor
        numerator = numerator.subs(F, sp.sympify(li.F))

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

        self.assertEqual( (sp.sympify(li.U) - target_U).simplify() , 0 )
        self.assertEqual( (sp.sympify(li.F) - target_F).simplify() , 0 )
        self.assertEqual( (sp.sympify(numerator) - target_numerator).simplify() , 0 )


        rank = 3
        target_exponent_U = N_nu - dim/2 * (L+1) - rank
        target_exponent_F = -(N_nu - dim/2 * L)
        target_exponentiated_F = target_F ** target_exponent_F
        target_exponentiated_U = target_U ** target_exponent_U
        self.assertEqual( (sp.sympify(li.exponentiated_U) - target_exponentiated_U).simplify() , 0 )
        self.assertEqual( (sp.sympify(li.exponentiated_F) - target_exponentiated_F).simplify() , 0 )

    #@attr('active')
    @attr('slow')
    def test_rank2_numerator_2L_box(self):
        k1, k2, p1, p2, p3, p4, s, t, mu, scalar_factor, D, eps = sp.symbols('k1 k2 p1 p2 p3 p4 s t mu scalar_factor D eps')

        z = sp.sympify(['z%i'%i for i in range(7+1)]) # Feynman parameters (only use z[1], z[2], ..., z[7] but not z[0])
        loop_momenta = [k1, k2]
        external_momenta = [p1, p2, p3, p4]
        indices = [mu]

        propagators = [k1**2,(k1+p2)**2,(k1-p1)**2,(k1-k2)**2,(k2+p2)**2,(k2-p1)**2,(k2+p2+p3)**2]
        numerator1 = 2*k1(mu)*k2(mu)
        numerator2 = 2*k1(mu)*k1(mu)
        replacement_rules = [(p1**2, 0),
                             (p2**2, 0),
                             (p3**2, 0),
                             (p4**2, 0),
                             (p1*p2, s/2),
                             (p2*p3, t/2),
                             (p1*p3, -s/2 - t/2)]

        li1 = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                            numerator=numerator1, Feynman_parameters=z[1:],
                                            Lorentz_indices=indices, dimensionality=D,
                                            replacement_rules=replacement_rules)
        li2 = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta,
                                            numerator=numerator2, Feynman_parameters=z[1:],
                                            Lorentz_indices=indices)
        Feynman_parametrized_numerator1 = sp.sympify(li1.numerator)
        Feynman_parametrized_numerator2 = sp.sympify(li2.numerator)

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

        self.assertEqual( (sp.sympify(li1.U) - target_U).simplify() , 0 )
        self.assertEqual( (sp.sympify(li2.U) - target_U).simplify() , 0 )
        self.assertEqual( (sp.sympify(li1.F) - target_F1).simplify() , 0 )
        self.assertEqual( (sp.sympify(li2.F) - target_F2).simplify() , 0 )

        for i, Feynman_parametrized_numerator in enumerate([Feynman_parametrized_numerator1, Feynman_parametrized_numerator2]):
            print(i)

            # need to insert the `F`s when comparing with SecDec
            Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs('F', target_Fs[i])

            # SecDec divides by a factor of ``Gamma(7 - D)`` while pySecDec divides by ``Gamma(6 - D)``
            Feynman_parametrized_numerator *= li1.Gamma_factor / sp.sympify('gamma(7 - D)')
            Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs('gamma(7 - D)', '(6 - D)*gamma(6 - D)')

            # SecDec puts the ``(-1)**N_nu`` into the prefactor while `LoopIntegral` puts it into the numerator.
            Feynman_parametrized_numerator = - Feynman_parametrized_numerator

            # only one scalar product per term --> can safely remove dummy index
            Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs(p1(mu), p1).subs(p2(mu), p2).subs(p3(mu), p3).subs(p4(mu), p4)

            if i == 1:
                Feynman_parametrized_numerator = Feynman_parametrized_numerator.subs(D, 4-2*eps)

            self.assertEqual( (sp.sympify(Feynman_parametrized_numerator) - target_numerators[i]).simplify() , 0 )

# TODO: uncomment when implemented
#    def test_error_if_index_too_often(self):
#        mu, nu = sp.symbols('mu nu')
#        numerator = 'k1(mu)*k2(mu)*k2(1)*k2(2)*p1(1)*p2(2)*p1(mu)'
#        loop_momenta = ['k1', 'k2']
#        external_momenta = ['p1', 'p2']
#        propagators = [] # dummy, do not need to specify propagators for the double index notation
#
#        li = LoopIntegral_from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator)
#        self.assertRaisesRegexp(AssertionError, '(E|e)ach.*index.*(exactly|at most).*(twice|two times)', lambda: li.numerator)
