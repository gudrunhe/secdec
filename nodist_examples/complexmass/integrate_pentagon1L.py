from __future__ import print_function
from pySecDec.integral_interface import IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load c++ library
    integral = IntegralLibrary('pentagon_1L/pentagon_1L_pylink.so')

    # choose integrator
    integral.use_Vegas(flags=2,epsrel=1e-4,epsabs=1e-13,nstart=10000,nincrease=1000,maxeval=100000000) # ``flags=2``: verbose --> see Cuba manual
    #integral.use_Divonne(flags=2,epsrel=1e-4,border=1e-7,maxeval=10000000) # ``flags=2``: verbose --> see Cuba manual

    number_of_real_parameters = int(integral.info['number_of_real_parameters'])
    number_of_complex_parameters = int(integral.info['number_of_complex_parameters'])

    with open('kinem.input') as f:
        with open('results_pentagon1L.txt', 'w') as resultsfile:
            for line in f:
                point = line.split()
                assert len(point) == 1+number_of_real_parameters+number_of_complex_parameters, "Invalid point: " + str(point)

                # convert to float and complex
                name = point[0]
                vals_real_parameters = [float(point[1+i]) for i in range(number_of_real_parameters)]
                vals_complex_parameters = [complex(point[1+number_of_real_parameters+i]) for i in range(number_of_complex_parameters)]

                # compute the integral
                str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral(vals_real_parameters,vals_complex_parameters)

                # print the result to resultsfile
                resultsfile.write('point: ' + str(point) + '\n')
                resultsfile.write('result: ' + str_integral_with_prefactor + '\n')
                resultsfile.write('--------\n\n\n')



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

    # examples how to access individual orders
    print('Numerical Result')
    print('prefactor=',prefactor)
    print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

    # results: checked with golem95

    # point: ['p1', '-9.025', '-2.', '-3.', '-4.', '-5.', '1.-1.5j', '1.2-0.2j', '2.5-0.3j']
    # result:  + ((7.75145915607855247e-03,4.71825754544246898e-03) +/- (1.94881592920598354e-06,1.74637246312473542e-06))*eps^-1 + ((-2.92107229713085836e-03,-3.66380586416755784e-03) +/- (4.36175447616102936e-06,4.03688104995995390e-06)) + O(eps)
    # --------
    # point: ['p2', '-9.025', '-2.', '-3.', '-4.', '-5.', '1.-2j', '1.2-0.2j', '2.5-0.3j']
    # result:  + ((6.74606643059131456e-03,5.09049971649342195e-03) +/- (1.94311745114479855e-06,1.69572076091542053e-06))*eps^-1 + ((-2.04653300788987941e-03,-3.47231822400529896e-03) +/- (4.26991608658848767e-06,3.89013102819715257e-06)) + O(eps)
    # --------
    # point: ['p3', '-9.025', '-2.', '-3.', '-4.', '-5.', '1.-2j', '1.2-2.2j', '2.5-3j']
    # result:  + ((4.15600555508943822e-04,6.44403609808399860e-03) +/- (1.44572946203320806e-06,1.17764875981478821e-06))*eps^-1 + ((1.29105597935993445e-03,-5.19095519363075783e-04) +/- (2.94003982584625173e-06,2.44375197454534635e-06)) + O(eps)
    # --------
