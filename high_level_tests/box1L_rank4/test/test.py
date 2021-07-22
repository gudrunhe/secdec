from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../box1L_rank4/box1L_rank4_pylink.so')
        self.real_parameters = [16.0, -75.0, 1.0]
        self.maxeval = 10**8
        self.epsrel = 1e-10
        self.epsabs = 1e-13

        self.target_result_without_prefactor = \
        {
              -1:   97.52083333333 +  0.0j,
               0: -181.22792152123 - 11.903058327787j
        }
        self.target_prefactor = \
        {
               0:  1.0
        }
        self.target_result_with_prefactor = \
        {
              -1:   97.52083333333 +  0.0j,
               0: -181.22792152123 - 11.903058327787j
        }

    def check_result(self, computed_series, target_series, epsrel, epsabs, order_min, order_max):
        # convert result to sympy expressions
        computed_series = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(order_min, order_max+1):
            value = complex( computed_series.coeff('eps',order).coeff('value') )
            error = complex( computed_series.coeff('eps',order).coeff('error') )

            # check that the uncertainties are reasonable
            self.assertLessEqual(error.real, abs(2*epsrel * target_series[order].real))
            if target_series[order].imag != 0.0:
                self.assertLessEqual(error.imag, abs(2*epsrel * target_series[order].imag))

            # check that the desired uncertainties are reached
            self.assertLessEqual(error.real, abs(epsrel * value.real) )
            if target_series[order].imag == 0.0:
                self.assertLessEqual(error.imag, epsabs)
            else:
                self.assertLessEqual(error.imag, abs(epsrel * value.imag) )

            # check integral value
            self.assertAlmostEqual(  value.real, target_series[order].real, delta=3.*epsrel*abs(target_series[order].real)  )
            if target_series[order].imag == 0.0:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=3.*epsabs  )
            else:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=3.*epsrel*abs(target_series[order].imag)  )

    def test_Cuhre(self):
        # choose integrator
        self.lib.use_Cuhre(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, real_complex_together=True)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters)

        # check integral
        #self.check_result(str_integral_without_prefactor, self.target_result_without_prefactor, self.epsrel, self.epsabs, order_min=0, order_max=1)
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs, order_min=-1, order_max=0)

        # check prefactor
        #prefactor = sp.sympify(  str_prefactor.replace(',','+I*')  )
        #for order in (-1,0):
        #    self.assertAlmostEqual(  prefactor.coeff('eps', order), self.target_prefactor[order], delta=1e-13  )

if __name__ == '__main__':
    unittest.main()
