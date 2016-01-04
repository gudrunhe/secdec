"""Unit tests for the U, F routines"""

from .loop_integral import *
import sympy as sp
import unittest
from nose.plugins.attrib import attr

def uf_from_propagators_generic(test_case, loop_momenta, propagators, result_u, result_f):
    loop_integral = LoopIntegral.from_propagators(propagators, loop_momenta, Feynman_parameters='Feynman')

    u,f = loop_integral.U, loop_integral.F

    sympy_u = sp.sympify(str(u))
    sympy_f = sp.sympify(str(f))

    result_u = sp.sympify(str(result_u).replace('x','Feynman'))
    result_f = sp.sympify(str(result_f).replace('x','Feynman'))

    zerou = (sympy_u - result_u).expand()
    zerof = (sympy_f - result_f).expand()

    test_case.assertEqual(zerou,0)
    test_case.assertEqual(zerof,0)

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
                                LoopIntegral.from_propagators, ['k1**2'], ['k1'], Feynman_parameters=['z0','z1'])

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
        loop_integral = LoopIntegral.from_propagators(propagators, loop_momenta, replacement_rules=replacement_rules, Feynman_parameters='z')
        U = sp.sympify(loop_integral.U)
        F = sp.sympify(loop_integral.F)

        target_U = sp.sympify('z3*(z4 + z5 + z6) + z0*(z3 + z4 + z5 + z6) + z1*(z3 + z4 + z5 + z6) + z2*(z3 + z4 + z5 + z6)')
        target_F = sp.sympify('-(s*z3*z4*z5) + s*z2*(-(z3*z4) - z4*z5) - z0*(s*z4*z5 + t*z3*z6) - s*z1*((z3 + z4)*z5 + z2*(z3 + z4 + z5 + z6))')

        self.assertEqual( (target_U - U).simplify() , 0)
        self.assertEqual( (target_F - F).simplify() , 0)

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
        loop_integral = LoopIntegral.from_propagators(props, loop_momenta, replacement_rules=replacement_rules, Feynman_parameters=['z1','z2','z3','z4'])
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
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegral.from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator)
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [
                            [(0,mu),(1,mu)], # (loop_momenta[0], index), (loop_momenta[1], index)
                            [(0,mu)],        # (loop_momenta[0], index)
                            [(1,mu),(1,nu)], # (loop_momenta[0], index), (loop_momenta[1], index)
                         ]

        self.assertEqual(len(tensors), 3)
        for i in range(3):
            print(i)
            self.assertEqual(tensors[i], target_tensors[i])

    #@attr('active')
    def test_double_index_notation_special_case(self):
        mu = sp.symbols('mu')
        numerator = 'k1(mu)'
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegral.from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator)
        tensors = li.numerator_loop_tensors

        target_tensors = [  [(0,mu)]  ]
        self.assertEqual(tensors, target_tensors)

    def test_double_index_notation_2(self):
        mu, nu = sp.symbols('mu nu')
        numerator = 'k1(mu)*k2(mu)*k2(1)*k2(2)*p1(1)*p2(2)'
        loop_momenta = ['k1', 'k2']
        external_momenta = ['p1', 'p2']
        propagators = [] # dummy, do not need to specify propagators for the double index notation

        li = LoopIntegral.from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator)
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [[(0,mu),(1,1),(1,2),(1,mu)]]

        self.assertEqual(len(tensors), 1)
        for i in range(1):
            print(i)
            self.assertEqual(tensors[i], target_tensors[i])

    #@attr('active')
    def test_tensor_numerator_bubble_1L(self): #TODO: test case for contracted loop momenta
        numerator = 'k(1)*k(2)*k(3)*k(4)'
        loop_momenta = ['k']
        external_momenta = ['p']
        propagators = ['k**2', '(k - p)**2']

        li = LoopIntegral.from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator, Feynman_parameters=['x1','x2'])
        tensors = li.numerator_loop_tensors
        # Note: `sympy` may reorder the terms
        target_tensors = [[(0,1),(0,2),(0,3),(0,4)]]

        self.assertEqual(len(tensors), 1)
        for i in range(1):
            print(i)
            self.assertEqual(tensors[i], target_tensors[i])

        numerator = li.numerator
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
                                      ''')
        self.assertEqual( (numerator - target_numerator).simplify() , 0 )

# TODO: uncomment when implemented
#    def test_error_if_index_too_often(self):
#        mu, nu = sp.symbols('mu nu')
#        numerator = 'k1(mu)*k2(mu)*k2(1)*k2(2)*p1(1)*p2(2)*p1(mu)'
#        loop_momenta = ['k1', 'k2']
#        external_momenta = ['p1', 'p2']
#        propagators = [] # dummy, do not need to specify propagators for the double index notation
#
#        li = LoopIntegral.from_propagators(propagators, loop_momenta, external_momenta, numerator=numerator)
#        self.assertRaisesRegexp(AssertionError, '(E|e)ach.*index.*(exactly|at most).*(twice|two times)', lambda: li.numerator)
