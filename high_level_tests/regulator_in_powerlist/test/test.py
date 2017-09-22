from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import unittest

class CheckLib(unittest.TestCase):
    def setUp(self):
        # load c++ library
        self.lib = IntegralLibrary('../regulator_in_powerlist/regulator_in_powerlist_pylink.so')

        # set global options
        self.maxeval = 10**7
        self.epsrel = 1e-3
        self.epsabs = 1e-4
        self.lib.use_Cuhre(self.epsrel, self.epsabs, flags=2,
                           real_complex_together=True,
                           maxeval=self.maxeval)
        self.lib.use_CQuad(self.epsrel, self.epsabs, verbose=True)

        # expected result
        self.target_result_without_prefactor_without_kinematics = sp.sympify('''
             + 7/3        * eps ** -4
             + 0          * eps ** -3
             - 30.7054359 * eps ** -2
        ''')
        self.target_prefactor = sp.sympify('''
                (-s)**(-4*eps - 1)*exp(-I*pi*(2*eps + 3))*gamma(4*eps + 1)/gamma(eps)**2
            ''').series(sp.symbols('eps'), n=5)
        self.order_min = -2
        self.order_max =  0

    def check_result(self, s):
        # compute integral
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = self.lib([s])

        # convert computed result to sympy expressions
        computed_series = sp.sympify(  str_integral_with_prefactor.replace(',','+I*').replace('+/-','*value+error*')  )

        # get target result (with nonstandard prescription of the negative log)
        log_minus_s = sp.log(-s) if s<0 else sp.log(s) - sp.I * sp.pi
        target_result_with_prefactor = \
            (
                self.target_result_without_prefactor_without_kinematics * \
                self.target_prefactor
            ).expand().replace('log(-s)',log_minus_s).subs('s',s)

        # print results
        print("obtained result without prefactor:\n", str_integral_without_prefactor, '\n')
        print("obtained result:\n", str_integral_with_prefactor, '\n')
        print("target result:\n", target_result_with_prefactor.evalf().expand(), '\n')
        print("----------------\n")

        for order in range(self.order_min, self.order_max+1):
            print('checking order "eps^%i"' % order)

            computed_value = complex( computed_series.coeff('eps',order).coeff('value') )
            computed_error = complex( computed_series.coeff('eps',order).coeff('error') )
            target_value = complex( target_result_with_prefactor.coeff('eps',order) )

            # check values (`2*epsrel` because the orders are mixed by the prefactor.)
            if target_value.real == 0.0:
                self.assertLessEqual(computed_error.real, self.epsabs)
                self.assertLessEqual(computed_value.real, 2. * computed_error.real)
            else:
                self.assertAlmostEqual(computed_value.real, target_value.real, delta=abs(2*self.epsrel*target_value.real))
            if target_value.imag == 0.0:
                self.assertLessEqual(computed_error.imag, self.epsabs)
                self.assertLessEqual(computed_value.imag, 2. * computed_error.imag)
            else:
                self.assertAlmostEqual(computed_value.imag, target_value.imag, delta=abs(2*self.epsrel*target_value.imag))

        print("----------------\n")

    def test_Euclidean_point(self):
        # check integral
        self.check_result(-2.0)

    def test_physical_point(self):
        # check integral
        self.check_result(+2.0)

if __name__ == '__main__':
    unittest.main()
