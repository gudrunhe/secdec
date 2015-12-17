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
