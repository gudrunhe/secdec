from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../userdefined_cpp/userdefined_cpp_pylink.so')
        self.target_result = {
                                 -2:   0.095206,
                                 -1:  -2.561,
                                  0:  21.120,
                                  1: -78.7
                             }
        self.real_parameters = [0.77]
        self.epsrel = 1e-3
        self.maxeval = 10**8
        self.epsabs_tol = 1e-15

    def check_result(self, computed_series, epsrel):
        # convert result to sympy expressions
        integral_with_prefactor = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(-2,2):
            value = complex( integral_with_prefactor.coeff('eps',order).coeff('value') )
            error = complex( integral_with_prefactor.coeff('eps',order).coeff('error') )

            # check that the uncertainties are reasonable
            self.assertLessEqual(error.real, abs(2*epsrel * self.target_result[order].real)+self.epsabs_tol)
            self.assertLessEqual(error.imag, abs(2*epsrel * self.target_result[order].imag)+self.epsabs_tol)

            # check that the desired uncertainties are reached
            self.assertLessEqual(error.real, abs(epsrel * value.real)+self.epsabs_tol )
            self.assertLessEqual(error.imag, abs(epsrel * value.imag)+self.epsabs_tol )

            # check integral value
            self.assertAlmostEqual(  value.real, self.target_result[order].real, delta=epsrel*abs(self.target_result[order].real)+self.epsabs_tol  )
            self.assertAlmostEqual(  value.imag, self.target_result[order].imag, delta=epsrel*abs(self.target_result[order].imag)+self.epsabs_tol  )

    def test_Vegas(self):
        # choose integrator
        self.lib.use_Vegas(flags=2, epsrel=self.epsrel, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters)

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Suave(self):
        # choose integrator
        self.lib.use_Suave(flags=2, epsrel=self.epsrel, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters)

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Divonne(self):
        # choose integrator
        self.lib.use_Divonne(flags=2, epsrel=self.epsrel, border=1e-8, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters)

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Cuhre(self):
        # choose integrator
        # Note: Need `mineval` because Cuhre underestimates the error
        self.lib.use_Cuhre(flags=2, epsrel=self.epsrel, maxeval=self.maxeval, mineval=6*10**4) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters)

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

if __name__ == '__main__':
    unittest.main()
