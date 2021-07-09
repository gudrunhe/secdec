from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../tadpole2L_rank2/tadpole2L_rank2_pylink.so')

        # set global options
        self.real_parameters = [1., 1., 1.]
        self.complex_parameters = []
        self.maxeval = 10**6
        self.epsrel = 1e-8
        self.epsabs = 1e-10

        self.target_result_without_prefactor = \
        {
              -1:  10.0,
               0: - 4.0,
               1: -34.5127841
        }

        self.target_prefactor = \
        {
              -1:  -0.25,
               0:  -0.461392167549234,
               1:  -1.87323249797567
        }
        
        self.target_result_with_prefactor = \
        {
              -2:  -2.5,
              -1:  -3.61392167549233578,
               0:  -8.25856028440502143
        }

    def check_result(self, computed_series, target_series, epsrel, epsabs, order_min, order_max):
        # convert result to sympy expressions
        if '+/-' in computed_series:
            computed_series = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )
        else:
            computed_series = (sp.sympify(  computed_series.replace(',','+I*')  ) * sp.sympify('value')).expand()

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
            self.assertAlmostEqual(  value.real, target_series[order].real, delta=epsrel*abs(target_series[order].real)  )
            if target_series[order].imag == 0.0:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=epsabs  )
            else:
                self.assertAlmostEqual(  value.imag, target_series[order].imag, delta=epsrel*abs(target_series[order].imag)  )

    def test_Cuhre(self):
        # choose integrator
        self.lib.use_Cuhre(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, real_complex_together=True)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters, self.complex_parameters)
        
        # check integral
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs, order_min=-2, order_max=0)

        # check prefactor
        #self.check_result(str_prefactor, self.target_prefactor, self.epsrel, self.epsabs, order_min=-2, order_max=0)

    def test_Cuhre_CQuad(self):
        # choose integrator
        self.lib.use_Cuhre(epsrel=self.epsrel, maxeval=self.maxeval, epsabs=self.epsabs, real_complex_together=True)
        self.lib.use_CQuad(epsrel=self.epsrel, epsabs=self.epsabs)

        # integrate
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib(self.real_parameters, self.complex_parameters)

        # check integral
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs, order_min=-2, order_max=0)

        # check prefactor
        #self.check_result(str_prefactor, self.target_prefactor, self.epsrel, self.epsabs, order_min=-2, order_max=0)

if __name__ == '__main__':
    unittest.main()
