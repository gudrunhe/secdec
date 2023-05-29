from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../difference/difference_pylink.so')
        self.target_result = {
                                  0:  1.64493406684822643647241516664602518923 + 3.1415926535897932384626433832795028842j,
                                  1:  2.08781123053685858754509217178101012328 - 6.2831853071795864769252867665590057684j,
                                  2: -5.94029019737039970544633397517750766917 + 4.2570651807194096861418776386549427857j
                             }
        self.epsrel = 1e-2

    def check_result(self, computed_series, epsrel):
        # convert result to sympy expressions
        integral_with_prefactor = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(3):
            value = complex( integral_with_prefactor.coeff('eps',order).coeff('value') )
            error = complex( integral_with_prefactor.coeff('eps',order).coeff('error') )

            # check that the uncertainties are reasonable
            self.assertLessEqual(error.real, abs(2*epsrel * self.target_result[order]))
            self.assertLessEqual(error.imag, abs(2*epsrel * self.target_result[order]))

            # check that the desired uncertainties are reached
            self.assertLessEqual(error.real, abs(epsrel * value) )
            self.assertLessEqual(error.imag, abs(epsrel * value) )

            # check integral value
            self.assertAlmostEqual(  value.real, self.target_result[order].real, delta=3.*epsrel*abs(self.target_result[order].real)  )
            self.assertAlmostEqual(  value.imag, self.target_result[order].imag, delta=3.*epsrel*abs(self.target_result[order].imag)  )

    def test_Vegas(self):
        # choose integrator
        self.lib.use_Vegas(flags=2, epsrel=self.epsrel) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Suave(self):
        # choose integrator
        self.lib.use_Suave(flags=2, epsrel=self.epsrel) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Divonne(self):
        # choose integrator
        self.lib.use_Divonne(flags=2, epsrel=self.epsrel, border=1e-8) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Cuhre(self):
        # choose integrator
        self.lib.use_Cuhre(flags=0, epsrel=self.epsrel) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

if __name__ == '__main__':
    unittest.main()
