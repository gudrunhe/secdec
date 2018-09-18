from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

# load c++ library
int5 = IntegralLibrary('crossint5/crossint5_pylink.so')

# choose integrator
#int5.use_Vegas(nstart=10**5, nincrease=5*10**4, epsrel=1e-3, epsabs=1e-7, maxeval=10**6, real_complex_together=True, flags=2) # ``flags=2``: verbose --> see Cuba manual
int5.use_Qmc(minn=10**8, minm=64, epsrel=1e-3, epsabs=1e-5, maxeval=10**10, verbosity=3, devices=[0,1,2,3], cudablocks=128, cudathreadsperblock=64, transform='korobov3')

number_of_real_parameters = int(int5.info['number_of_real_parameters'])
number_of_complex_parameters = int(int5.info['number_of_complex_parameters'])


with open('kinematics.input') as f:
    with open('results_crossint5.out', 'w') as resultsfile:
        for line in f:
            point = line.split()
            assert len(point) == 1+number_of_real_parameters+number_of_complex_parameters, "Invalid point: " + str(point)

            # convert to float and complex
            name = point[0]
            vals_real_parameters = [float(point[1+i]) for i in range(number_of_real_parameters)]
            vals_complex_parameters = [complex(point[1+number_of_real_parameters+i]) for i in range(number_of_complex_parameters)]

            # compute the integral
            str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = \
                int5(
                        vals_real_parameters,
                        vals_complex_parameters,
                        number_of_presamples=10**5,
                        deformation_parameters_maximum=1e-2,
                        together=True
                    )

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
