from . import integral_interface as ii
import unittest

class TestParsing(unittest.TestCase):
    def test_parse_series_coefficient(self):
        self.assertEqual(ii._parse_series_coefficient("1"), 1)
        self.assertEqual(ii._parse_series_coefficient("1 +/- 2"), (1, 2))
        self.assertEqual(ii._parse_series_coefficient("(1,2)"), 1 + 2j)
        self.assertEqual(ii._parse_series_coefficient("(1,2) +/- (3,4)"), (1 + 2j, 3 + 4j))

    def test_parse_series(self):
        self.assertEqual(
            ii._parse_series(" + O(ep^2)"),
            ([], "ep", 2)
        )
        self.assertEqual(
            ii._parse_series(" + (1) + O(ep^2)"),
            ([(1, 0)], "ep", 2)
        )
        self.assertEqual(
            ii._parse_series(" + (1.25000000e+00)*ep + O(ep^2)"),
            ([(1.25, 1)], "ep", 2)
        )
        self.assertEqual(
            ii._parse_series(" + ((1.25000000e+00,2.5e+01)) + O(ep)"),
            ([(1.25+25j, 0)], "ep", 1)
        )
        self.assertEqual(
            ii._parse_series(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)"),
            ([((0.5+1j, 2+3j), -1), ((4+5j, 6+7j), 1)], "ep", 2)
        )
        self.assertEqual(
            ii._parse_series(" + (1) + (-2)*ep"),
            ([(1,0),(-2,1)], "ep", None)
        )

    def test_series_to_ginac(self):
        self.assertEqual(
            ii.series_to_ginac(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)"),
            ("(0.5+1*I)/ep + (4+5*I)*ep + Order(ep^2)", "(2+3*I)/ep + (6+7*I)*ep + Order(ep^2)")
        )
        self.assertEqual(
            ii.series_to_ginac(" + (1)*ep^-1 + (2)"),
            ("1/ep + 2", "0/ep + 0")
        )
        self.assertEqual(
            ii.series_to_ginac(" + (8)"),
            ("8", "0")
        )

    def test_series_to_sympy(self):
        self.assertEqual(
            ii.series_to_sympy(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)"),
            ("(0.5+1*I)/ep + (4+5*I)*ep + O(ep**2)", "(2+3*I)/ep + (6+7*I)*ep + O(ep**2)")
        )
        self.assertEqual(
            ii.series_to_sympy(" + (1)*ep^-1 + (2)"),
            ("1/ep + 2", "0/ep + 0")
        )
        self.assertEqual(
            ii.series_to_sympy(" + (8)"),
            ("8", "0")
        )

    def test_series_to_mathematica(self):
        self.assertEqual(
            ii.series_to_mathematica(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)"),
            ("(0.5+1*I)/ep + (4+5*I)*ep + O[ep]^2", "(2+3*I)/ep + (6+7*I)*ep + O[ep]^2")
        )
        self.assertEqual(
            ii.series_to_mathematica(" + (1)*ep^-1 + (2)"),
            ("1/ep + 2", "0/ep + 0")
        )
        self.assertEqual(
            ii.series_to_mathematica(" + (8)"),
            ("8", "0")
        )

    def test_series_to_maple(self):
        self.assertEqual(
            ii.series_to_maple(" + ((5e-01,10e-01) +/- (20e-01,3e-00))*ep^-1 + ((4,5) +/- (6,7))*ep + O(ep^2)"),
            ("(0.5+1*I)/ep + (4+5*I)*ep + O(ep^2)", "(2+3*I)/ep + (6+7*I)*ep + O(ep^2)")
        )
        self.assertEqual(
            ii.series_to_maple(" + (1)*ep^-1 + (2)"),
            ("1/ep + 2", "0/ep + 0")
        )
        self.assertEqual(
            ii.series_to_maple(" + (8)"),
            ("8", "0")
        )
