from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
triangle = IntegralLibrary('triangle2L/triangle2L_pylink.so')

# choose integrator
triangle.use_Vegas(flags=2) # ``flags=2``: verbose --> see Cuba manual

number_of_real_parameters = int(triangle.info['number_of_real_parameters'])
number_of_complex_parameters = int(triangle.info['number_of_complex_parameters'])

with open('kinematics.input') as f:
  with open('results_triangle2L.txt', 'w') as resultsfile:
    for line in f:
        point = line.split()
        assert len(point) == 1+number_of_real_parameters+number_of_complex_parameters, "Invalid point: " + str(point)

        # convert to float and complex
        name = point[0]
        vals_real_parameters = [float(point[1+i]) for i in range(number_of_real_parameters)]
        vals_complex_parameters = [complex(point[1+number_of_real_parameters+i]) for i in range(number_of_complex_parameters)]

        # compute the integral
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = triangle(vals_real_parameters,vals_complex_parameters)

        # print numerical point read from file kinematics.input:
        print('read point "' + name + '"')
        print('real_parameters:', vals_real_parameters)
        print('complex_parameters:', vals_complex_parameters)

        # convert complex numbers from c++ to sympy notation
        str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
        str_prefactor = str_prefactor.replace(',','+I*')
        str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')

        # convert result to sympy expressions
        integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
        integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
        prefactor = sp.sympify(str_prefactor)
        integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
        integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

        # example how to access individual orders
        print('leading pole:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
        print('subleading pole:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
        print('finite part:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

        # example how to print the result to a file
        resultsfile.write('point: ' + str(point) + '\n')
        resultsfile.write('leading_pole = ' + str(integral_with_prefactor.coeff('eps',-2).coeff('value')) + '\n')
        resultsfile.write('err_leading_pole = ' + str(integral_with_prefactor_err.coeff('eps',-2).coeff('error')) + '\n')
        resultsfile.write('subleading_pole = ' + str(integral_with_prefactor.coeff('eps',-1).coeff('value')) + '\n')
        resultsfile.write('err_subleading_pole = ' + str(integral_with_prefactor_err.coeff('eps',-1).coeff('error')) + '\n')
        resultsfile.write('finite_part = ' + str(integral_with_prefactor.coeff('eps',0).coeff('value')) + '\n')
        resultsfile.write('err_finite_part = ' + str(integral_with_prefactor_err.coeff('eps',0).coeff('error')) + '\n')
        resultsfile.write('--------\n\n\n')
