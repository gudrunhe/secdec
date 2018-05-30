from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckQmcErrorMessages(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../one_integration_variable/one_integration_variable_pylink.so')

    def test_setting_errormode(self):
        self.assertRaisesRegexp(ValueError, 'Unknown `errormode` "foo"', self.lib.use_Qmc, errormode='foo')

        # test known errormodes
        self.lib.use_Qmc(errormode='default')
        self.lib.use_Qmc(errormode='all')
        self.lib.use_Qmc(errormode='largest')

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../one_integration_variable/one_integration_variable_pylink.so')
        self.target_result = {
                                 -1:  0.0,
                                  0: -0.5,
                                  1: -0.25
                             }
        self.epsrel = 1e-11
        self.epsabs = 1e-10
        self.maxeval = 10**5

    def check_result(self, computed_series, epsrel, epsabs):
        # convert result to sympy expressions
        integral_with_prefactor = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(-1,2):
            value = complex( integral_with_prefactor.coeff('eps',order).coeff('value') )
            error = complex( integral_with_prefactor.coeff('eps',order).coeff('error') )

            # check that the desired uncertainties are reached
            self.assertLessEqual(error.real, epsabs)
            self.assertLessEqual(error.imag, epsabs)

            # check integral value
            self.assertAlmostEqual(  value.real, self.target_result[order].real, delta=epsabs  )
            self.assertAlmostEqual(  value.imag, self.target_result[order].imag, delta=epsabs  )

    def test_Vegas(self):
        # choose integrator
        # can only reach ~2e-9 accuracy with Vegas
        self.lib.use_Vegas(flags=0, epsrel=self.epsrel, epsabs=2e-9, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel, 2e-9)

    def test_Suave(self):
        # choose integrator
        self.lib.use_Suave(flags=0, epsrel=self.epsrel, epsabs=self.epsabs, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel, self.epsabs)

    def test_Divonne(self):
        # choose integrator
        self.lib.use_Divonne(flags=0, epsrel=self.epsrel, epsabs=self.epsabs, border=1e-8, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel, self.epsabs)

    def test_Cuhre(self):
        # choose integrator
        self.lib.use_Cuhre(flags=0, epsrel=self.epsrel, epsabs=self.epsabs, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel, self.epsabs)

    def test_CQuad(self):
        # choose integrator
        self.lib.use_CQuad(verbose=False, epsrel=self.epsrel, epsabs=self.epsabs)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel, self.epsabs)

    def test_Qmc(self):
        # choose integrator
        self.lib.use_Qmc(verbosity=0, epsrel=self.epsrel, epsabs=self.epsabs, seed=3212)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel, self.epsabs)

if __name__ == '__main__':
    unittest.main()
