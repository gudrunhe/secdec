from pySecDec.integral_interface import IntegralLibrary
from math import log

# load c++ library
easy = IntegralLibrary('easy/easy_pylink.so')

# integrate
_, _, result = easy()

# print result
print('Numerical Result:' + result)
print('Analytic Result:' + ' + (%f)*eps^-1 + (%f) + O(eps)' % (1.0,1.0-log(2.0)))
