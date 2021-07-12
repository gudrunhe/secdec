from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../userdefined_cpp/userdefined_cpp_pylink.so')

        # full analytic result for func(y)=1:
        #   (4^-eps (-2 + 2^eps))/((-1 + eps) eps^2)
        # = 1/eps^2 + (1 - Log[2] - Log[4])/eps +  1/2 (2 - 2 Log[2] - Log[2]^2 - 2 Log[4] + 2 Log[2] Log[4] + Log[4]^2) + O[eps]
        # = 1.0000000000000000000/eps^2 - 1.0794415416798359283/eps + 0.60214400703386905808 + O[eps]
        self.target_result = {
                                 -2:   1.0,
                                 -1:  -1.0794415416798359283,
                                  0:   0.60214400703386905808
                             }
        self.epsrel = 1e-4
        self.maxeval = 10**8
        self.epsabs_tol = 1e-15

    def check_result(self, computed_series, epsrel):
        # convert result to sympy expressions
        integral_with_prefactor = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(-2,0):
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
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Suave(self):
        # choose integrator
        self.lib.use_Suave(flags=2, epsrel=self.epsrel, maxeval=self.maxeval, nnew=10**5) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Divonne(self):
        # choose integrator
        self.lib.use_Divonne(flags=2, epsrel=self.epsrel, border=1e-12, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

    def test_Cuhre(self):
        # choose integrator
        self.lib.use_Cuhre(flags=2, epsrel=self.epsrel, maxeval=self.maxeval) # ``flags=2``: verbose --> see Cuba manual

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib()

        # check
        self.check_result(str_integral_with_prefactor, self.epsrel)

if __name__ == '__main__':
    unittest.main()
