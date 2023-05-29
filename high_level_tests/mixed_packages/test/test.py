from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../mixed_packages/mixed_packages_pylink.so')

        # set global options
        self.real_parameters = [10., 2.]
        self.complex_parameters = []
        self.maxeval = 10**7
        self.epsrel = 1e-4
        self.epsabs = 1e-7
        
        # s/msq * bub(s,msq,msq) + msq/s* Integrate[1/(s z1 + msq z2 z3), {z1, 0, 1}, {z2, 0, 1}, {z3, 0, 1}]
        self.target_result_with_prefactor = \
        {
              -1:  5.0              +  0.0j,
               0: 10.63956085778167 +  7.024814731040726j,
        }
        self.order_min = -1
        self.order_max =  0

    def check_result(self, computed_series, target_series, epsrel, epsabs):
        # convert result to sympy expressions
        computed_series = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(self.order_min, self.order_max+1):
            value = complex( computed_series.coeff('eps',order).coeff('value') )
            error = complex( computed_series.coeff('eps',order).coeff('error') )

            # check that the uncertainties are reasonable
            self.assertLessEqual(error.real, abs(2*epsrel * target_series[order].real))
            if target_series[order].imag != 0.0:
                self.assertLessEqual(error.imag, abs(2*epsrel * target_series[order].imag))

            # check that the desired uncertainties are reached
            self.assertLessEqual(error.real, abs(epsrel * value.real) )
            if target_series[order].imag == 0.0:
                self.assertLessEqual(error.imag, 5*epsabs)
            else:
                self.assertLessEqual(error.imag, abs(epsrel * 5 * value.imag) )

            # check integral value
            self.assertAlmostEqual(  value.real, target_series[order].real, delta=3.*epsrel*abs(target_series[order].real)  )
            if target_series[order].imag == 0.0:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=3.*epsabs  )
            else:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=3.*epsrel*abs(target_series[order].imag)  )

    def test_Cuhre(self):
        # choose integrator
        self.lib.use_Cuhre(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, real_complex_together=True, flags=0)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters, self.complex_parameters)

        # check integral
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs)

    def test_Cuhre_CQuad(self):
        # choose integrator
        self.lib.use_Cuhre(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, real_complex_together=True, flags=0)
        self.lib.use_CQuad(epsrel=self.epsrel, epsabs=self.epsabs, verbose=False)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters, self.complex_parameters)

        # check integral
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs)

    def test_Qmc_default_integral_transform(self):
        # choose integrator
        self.lib.use_Qmc(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, verbosity=0, seed=143, transform='korobov3')

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters, self.complex_parameters)

        # check integral
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs)

if __name__ == '__main__':
    unittest.main()
