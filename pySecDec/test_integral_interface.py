from . import integral_interface as ii
import unittest
import pytest

#@pytest.mark.active
class TestParsing(unittest.TestCase):
    def test_parse_series_coefficient(self):
        assert ii._parse_series_coefficient("1") == 1
        assert ii._parse_series_coefficient("1 +/- 2") == (1, 2)
        assert ii._parse_series_coefficient("(1,2)") == 1 + 2j
        assert ii._parse_series_coefficient("(1,2) +/- (3,4)") == (1 + 2j, 3 + 4j)

    def test_parse_series(self):
        assert ii._parse_series(" + O(ep^2)") == \
            ([], "ep", 2)
        assert ii._parse_series(" + (1) + O(ep^2)") == \
            ([(1, 0)], "ep", 2)
        assert ii._parse_series(" + (1.25000000e+00)*ep + O(ep^2)") == \
            ([(1.25, 1)], "ep", 2)
        assert ii._parse_series(" + ((1.25000000e+00,2.5e+01)) + O(ep)") == \
            ([(1.25+25j, 0)], "ep", 1)
        assert ii._parse_series(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)") == \
            ([((0.5+1j, 2+3j), -1), ((4+5j, 6+7j), 1)], "ep", 2)
        assert ii._parse_series(" + (1) + (-2)*ep") == \
            ([(1,0),(-2,1)], "ep", None)

    def test_series_to_ginac(self):
        assert ii.series_to_ginac(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)") == \
            ("(0.5+1*I)/ep + (4+5*I)*ep + Order(ep^2)", "(2+3*I)/ep + (6+7*I)*ep + Order(ep^2)")
        assert ii.series_to_ginac(" + (1)*ep^-1 + (2)") == \
            ("1/ep + 2", "0/ep + 0")
        assert ii.series_to_ginac(" + (8)") == \
            ("8", "0")
        assert ii.series_to_ginac(" + ( + (-1 +/- 2)*eps^-1 + (-10 +/- 8) + O(eps))*alpha^-1 + ( + (5 +/- 1)*eps^-2 + (-19 +/- 9) + O(eps)) + O(alpha)") == \
            (
                "(-1/eps + -10 + Order(eps))/alpha + (5/eps^2 + -19 + Order(eps)) + Order(alpha)",
                "(2/eps + 8 + Order(eps))/alpha + (1/eps^2 + 9 + Order(eps)) + Order(alpha)"
            )

    def test_series_to_sympy(self):
        assert ii.series_to_sympy(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)") == \
            ("(0.5+1*I)/ep + (4+5*I)*ep + O(ep**2)", "(2+3*I)/ep + (6+7*I)*ep + O(ep**2)")
        assert ii.series_to_sympy(" + (1)*ep^-1 + (2)") == \
            ("1/ep + 2", "0/ep + 0")
        assert ii.series_to_sympy(" + (8)") == \
            ("8", "0")
        assert ii.series_to_sympy(" + ( + (-1 +/- 2)*eps^-1 + (-10 +/- 8) + O(eps))*alpha^-1 + ( + (5 +/- 1)*eps^-2 + (-19 +/- 9) + O(eps)) + O(alpha)") == \
            (
                "(-1/eps + -10 + O(eps))/alpha + (5/eps**2 + -19 + O(eps)) + O(alpha)",
                "(2/eps + 8 + O(eps))/alpha + (1/eps**2 + 9 + O(eps)) + O(alpha)"
            )

    def test_series_to_mathematica(self):
        assert ii.series_to_mathematica(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)") == \
            ("(0.5+1*I)/ep + (4+5*I)*ep + O[ep]^2", "(2+3*I)/ep + (6+7*I)*ep + O[ep]^2")
        assert ii.series_to_mathematica(" + (1)*ep^-1 + (2)") == \
            ("1/ep + 2", "0/ep + 0")
        assert ii.series_to_mathematica(" + (8)") == \
            ("8", "0")
        assert ii.series_to_mathematica(" + (1e-10 +/- 2e+20)") == \
            ("1.00000000000000004*10^-10", "2*10^+20")
        assert ii.series_to_mathematica(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)") == \
            ("(0.5+1*I)/ep + (4+5*I)*ep + O[ep]^2", "(2+3*I)/ep + (6+7*I)*ep + O[ep]^2")
        assert ii.series_to_mathematica(" + ( + (-1 +/- 2)*eps^-1 + (-10 +/- 8) + O(eps))*alpha^-1 + ( + (5 +/- 1)*eps^-2 + (-19 +/- 9) + O(eps)) + O(alpha)") == \
            (
                "SeriesData[alpha, 0, {-1/eps + -10 + O[eps], 5/eps^2 + -19 + O[eps]}, -1, 1, 1]",
                "SeriesData[alpha, 0, {2/eps + 8 + O[eps], 1/eps^2 + 9 + O[eps]}, -1, 1, 1]"
            )

    def test_series_to_maple(self):
        assert ii.series_to_maple(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)") == \
            ("(0.5+1*I)/ep + (4+5*I)*ep + O(ep^2)", "(2+3*I)/ep + (6+7*I)*ep + O(ep^2)")
        assert ii.series_to_maple(" + (1)*ep^-1 + (2)") == \
            ("1/ep + 2", "0/ep + 0")
        assert ii.series_to_maple(" + (8)") == \
            ("8", "0")
        assert ii.series_to_maple(" + ( + (-1 +/- 2)*eps^-1 + (-10 +/- 8) + O(eps))*alpha^-1 + ( + (5 +/- 1)*eps^-2 + (-19 +/- 9) + O(eps)) + O(alpha)") == \
            (
                "(-1/eps + -10 + O(eps))/alpha + (5/eps^2 + -19 + O(eps)) + O(alpha)",
                "(2/eps + 8 + O(eps))/alpha + (1/eps^2 + 9 + O(eps)) + O(alpha)"
            )

    def test_multi_series_to_sympy(self):
        assert ii.series_to_sympy(" + (1)*ep^-1 + (2)\n + (1)*ep^-2 + (2)*ep^-1") == \
            [('1/ep + 2', '0/ep + 0'), ('1/ep**2 + 2/ep', '0/ep**2 + 0/ep')]

    def test_multi_series_to_mathematica(self):
        assert ii.series_to_mathematica(" + (1)*ep^-1 + (2)\n + (1)*ep^-2 + (2)*ep^-1") == \
            [('1/ep + 2', '0/ep + 0'), ('1/ep^2 + 2/ep', '0/ep^2 + 0/ep')]

    def test_disteval_series_to_x(self):
        series = """[
          (
            +eta^0*eps^0*(+1.0e+00+0.0e+00j)
            +eta^0*eps^0*(+2.0e+00+0.0e+00j)*plusminus
            +eta^0*eps^1*(-3.0e+00+0.0e+00j)
            +eta^0*eps^1*(+4.0e+00+0.0e+00j)*plusminus
            +eta^1*eps^0*(+5.0e+00-9.0e+00j)
            +eta^1*eps^0*(+6.0e+00+0.0e+00j)*plusminus
            +eta^1*eps^1*(-7.0e+00-1.0e+01j)
            +eta^1*eps^1*(+8.0e+00+0.0e+00j)*plusminus
          )
        ]
        """
        assert ii.series_to_sympy(series) == \
            ('(1 + -3*eps + O(eps**2)) + ((5-9*I) + (-7-10*I)*eps + O(eps**2))*eta + O(eta**2)',
             '(2 + 4*eps + O(eps**2)) + (6 + 8*eps + O(eps**2))*eta + O(eta**2)')
        assert ii.series_to_mathematica(series) == \
            ('SeriesData[eta, 0, {1 + -3*eps + O[eps]^2, (5-9*I) + (-7-10*I)*eps + O[eps]^2}, 0, 2, 1]',
             'SeriesData[eta, 0, {2 + 4*eps + O[eps]^2, 6 + 8*eps + O[eps]^2}, 0, 2, 1]')

    def test_disteval_multi_series_to_x(self):
        series = """[
          (
            +eps^0*(+1.0e+01+2.0e+01j)
            +eps^0*(+3.0e+00+4.0e+00j)*plusminus
            +eps^1*(+5.0e+01+6.0e+01j)
            +eps^1*(+7.0e+00+8.0e+00j)*plusminus
          ),
          (
            +eps^1*(-9.0e+01-1.0e+01j)
            +eps^1*(-2.0e+00-3.0e+00j)*plusminus
          )
        ]
        """
        assert ii.series_to_sympy(series) == \
            [('(10+20*I) + (50+60*I)*eps + O(eps**2)',
              '(3+4*I) + (7+8*I)*eps + O(eps**2)'),
             ('(-90-10*I)*eps + O(eps**2)',
              '(-2-3*I)*eps + O(eps**2)')]
        assert ii.series_to_mathematica(series) == \
            [('(10+20*I) + (50+60*I)*eps + O[eps]^2',
              '(3+4*I) + (7+8*I)*eps + O[eps]^2'),
             ('(-90-10*I)*eps + O[eps]^2',
              '(-2-3*I)*eps + O[eps]^2')]
