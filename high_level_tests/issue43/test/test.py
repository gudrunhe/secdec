from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary, DistevalLibrary
import sympy as sp
import sys
import subprocess
import json
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ libraries
        self.lib = IntegralLibrary('../issue43/issue43_pylink.so')
        self.distlib = DistevalLibrary('../issue43/disteval/issue43.json')

        # set global options
        self.real_parameters = []
        self.complex_parameters = []
        self.maxeval = 10**6
        self.epsrel = 1e-2
        self.epsabs = 1e-5

        self.target_result_with_prefactor = \
        {
               -1: -1.0 + 0.0j,
               0: 0 + 0.0j,
        }
        self.order_min = -1
        self.order_max = 0

    def check_result(self, computed_series, target_series, epsrel, epsabs, order_min ,order_max, integrator='intlib'):
    
        if integrator == 'intlib':
            # convert result to sympy expressions
            computed_series = sp.sympify(  computed_series.replace(',','+I*').replace('+/-','*value+error*')  )

        for i, order in enumerate(range(order_min, order_max+1)):
        
            if integrator == 'intlib':
                value = complex( computed_series.coeff('eps',order).coeff('value') )
                error = complex( computed_series.coeff('eps',order).coeff('error') )
            elif integrator == 'disteval':
                value = computed_series[(order,)][0]
                error = computed_series[(order,)][1]
            elif integrator == 'disteval-cli':
                self.assertEqual(computed_series[i][0][0],order)
                value = complex(computed_series[i][1][0],computed_series[i][1][1])
                error = complex(computed_series[i][2][0],computed_series[i][2][1])
            else:
                raise ValueError(f"unknown integrator in check_result, got {integrator}")

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
        
        # check integral
        self.check_result(str_integral_with_prefactor, self.target_result_with_prefactor, self.epsrel, self.epsabs, self.order_min, self.order_max)
        
    def test_disteval(self):
        # integrate
        result = self.distlib(parameters={}, epsrel=self.epsrel, epsabs=self.epsabs, format='json', verbose=False)

        # check integral
        print(result)

        self.check_result(result["sums"]["issue43"], self.target_result_with_prefactor, self.epsrel, self.epsabs, self.order_min, self.order_max, integrator='disteval')
        
        
    def test_disteval_cli(self):
        psd_out = subprocess.run([f'{sys.executable} -m pySecDec.disteval ../issue43/disteval/issue43.json --epsabs={self.epsabs} --epsrel={self.epsrel} --format="json"'], shell=True, check=True, capture_output=True)
        result = json.loads(psd_out.stdout)
        
        # check integral
        self.check_result(result["sums"]["issue43"], self.target_result_with_prefactor, self.epsrel, self.epsabs, self.order_min, self.order_max, integrator='disteval-cli')

if __name__ == '__main__':
    unittest.main()
