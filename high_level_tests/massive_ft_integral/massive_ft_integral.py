#!/usr/bin/env python3
from pySecDec import sum_package, MakePackage

if __name__ == "__main__":

    poly_lower = ['(x1 + 1)**(eps - 2)']
    remainder_expression_lower = 'exp(-m*x1)'

    poly_upper = ['x1**(-eps)', '(x1 + 1)**(eps - 2)']
    remainder_expression_upper = 'exp(-m/x1)'

    lower_integral = MakePackage('lower_integral', real_parameters = ['t', 'm'], regulators = ['eps'], requested_orders = [3], integration_variables = ['x1'],
                polynomials_to_decompose = poly_lower, remainder_expression = remainder_expression_lower)
    upper_integral = MakePackage('upper_integral', real_parameters = ['t', 'm'], regulators = ['eps'], requested_orders = [3], integration_variables = ['x1'],
                polynomials_to_decompose = poly_upper, remainder_expression = remainder_expression_upper)

    sum_package('massive_ft_integral', [lower_integral, upper_integral], requested_orders = [3], real_parameters = ['t', 'm'], regulators = ['eps'])