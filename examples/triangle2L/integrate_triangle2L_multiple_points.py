from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
triangle2L = IntegralLibrary('triangle2L/triangle2L_pylink.so')

# choose integrator
triangle2L.use_Vegas(flags=2) # ``flags=2``: verbose --> see Cuba manual

number_of_real_parameters = int(triangle2L.info['number_of_real_parameters'])
number_of_complex_parameters = int(triangle2L.info['number_of_complex_parameters'])

with open('kinematics.input') as f:
  with open('results_P126.txt', 'w') as resultsfile:
    for line in f:
        point = line.split()
        assert len(point) == 1+number_of_real_parameters+number_of_complex_parameters, "Invalid point: " + str(point)

        # convert to float and complex
        name = point[0]
        vals_real_parameters = [float(point[1+i]) for i in range(number_of_real_parameters)]
        vals_complex_parameters = [complex(point[1+number_of_real_parameters+i]) for i in range(number_of_complex_parameters)]

        # compute the integral
        str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = triangle2L(vals_real_parameters,vals_complex_parameters)


        # print the result to resultsfile
        resultsfile.write('point: ' + str(point) + '\n')
        resultsfile.write('result: ' + str_integral_with_prefactor + '\n')
        resultsfile.write('--------\n\n\n')

        # alternative output format
        #resultsfile.write('point: ' + str(point) + '\n')
        #resultsfile.write('leading_pole = ' + str(integral_with_prefactor.coeff('eps',-2).coeff('value')) + '\n')
        #resultsfile.write('err_leading_pole = ' + str(integral_with_prefactor_err.coeff('eps',-2).coeff('error')) + '\n')
        #resultsfile.write('subleading_pole = ' + str(integral_with_prefactor.coeff('eps',-1).coeff('value')) + '\n')
        #resultsfile.write('err_subleading_pole = ' + str(integral_with_prefactor_err.coeff('eps',-1).coeff('error')) + '\n')
        #resultsfile.write('finite_part = ' + str(integral_with_prefactor.coeff('eps',0).coeff('value')) + '\n')
        #resultsfile.write('err_finite_part = ' + str(integral_with_prefactor_err.coeff('eps',0).coeff('error')) + '\n')
        #resultsfile.write('--------\n\n\n')
