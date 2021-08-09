from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../easy_sum/easy_sum_pylink.so')

        # set global options
        self.real_parameters = [0.1]
        self.complex_parameters = []
        self.maxeval = 10**6
        self.epsrel = 1e-4
        self.epsabs = 1e-7

        self.target_result_with_prefactor_1 = \
        {
              -1:                     0.2 + 0.0j, # eps ** -1
               0:  0.22962348064032504712 + 0.0j, # eps **  0
        }
        self.order_min_1 =  -1
        self.order_max_1 =  0
        
        self.target_result_with_prefactor_2 = \
        {
              -2:                     0.05 + 0.0j, # eps ** -2
              -1:  0.015342640972002734529 + 0.0j, # eps ** -1
               0: 0.0033313156240476989125 + 0.0j, # eps **  0
        }
        self.order_min_2 =  -2
        self.order_max_2 =  0

    def check_result(self, computed_series, target_series, epsrel, epsabs, order_min ,order_max):
        # convert result to sympy expressions
        computed_series = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for order in range(order_min, order_max+1):
            value = complex( computed_series.coeff('eps',order).coeff('value') )
            error = complex( computed_series.coeff('eps',order).coeff('error') )

            # check that the uncertainties are reasonable
            if target_series[order].real != 0.0:
                self.assertLessEqual(error.real, abs(2*epsrel * target_series[order].real))
            if target_series[order].imag != 0.0:
                self.assertLessEqual(error.imag, abs(2*epsrel * target_series[order].imag))

            # check that the desired uncertainties are reached
            if target_series[order].real == 0.0:
                self.assertLessEqual(error.real, max( abs(epsrel * value), epsabs) )
            else:
                self.assertLessEqual(error.real, abs(epsrel * value.real) )
            if target_series[order].imag == 0.0:
                self.assertLessEqual(error.imag, max( abs(epsrel * value), epsabs) )
            else:
                self.assertLessEqual(error.imag, abs(epsrel * value.imag) )

            # check integral value
            if target_series[order].real == 0.0:
                self.assertAlmostEqual(  value.real, target_series[order].real, delta=max( 3.*epsrel*abs(target_series[order]), 3.*epsabs)  )
            else:
                self.assertAlmostEqual(  value.real, target_series[order].real, delta=3.*epsrel*abs(target_series[order].real)  )
            if target_series[order].imag == 0.0:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=max( 3.*epsrel*abs(target_series[order]), 3.*epsabs)  )
            else:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=3.*epsrel*abs(target_series[order].imag)  )

    def test_Qmc(self):
        # choose integrator
        self.lib.use_Qmc(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, verbosity=0, seed=143, transform='korobov3')

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters, self.complex_parameters)
                
        str_integral_with_prefactor_1, str_integral_with_prefactor_2 = str_integral_with_prefactor.strip().split('\n')
        
        # check integral
        self.check_result(str_integral_with_prefactor_1, self.target_result_with_prefactor_1, self.epsrel, self.epsabs, self.order_min_1, self.order_max_1)
        self.check_result(str_integral_with_prefactor_2, self.target_result_with_prefactor_2, self.epsrel, self.epsabs, self.order_min_2, self.order_max_2)


if __name__ == '__main__':
    unittest.main()
