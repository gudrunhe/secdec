from pySecDec.integral_interface import IntegralLibrary
from math import log

# load c++ library
easy = IntegralLibrary('easy/easy_pylink.so')

# choose Qmc integrator
# automatically uses all avaliable GPUs
easy.use_Qmc(transform='korobov3')

# integrate
_, _, result = easy()

# print result
print('Numerical Result:' + result)
print('Analytic Result:' + ' + (%.15g)*eps^-1 + (%.15g) + O(eps)' % (1.0,1.0-log(2.0)))
