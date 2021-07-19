#!/usr/bin/env python3

from pySecDec.integral_interface import IntegralLibrary
import sympy as sp
import numpy as np

# define analytic result
def f(s,t,u, mtsq):
    """
    The resulting function of the integral from arXiv:1812.04373, to order 0 in mtsq.

    Expected Result from the paper for ['s','t','u','mtsq'] = [4.0, -2.82842712475, -2.82842712475, 0.1]:
    (-1.30718609739+1.85618421207j)
    """
    v = -t/s
    lsm = np.log(complex(s/mtsq))
    return 1/(s**2*v)*(np.pi**2-2*(lsm-1j*np.pi)*(lsm+np.log(complex(v))))

# define useful functions
def fourvec(m,p):
    # creates a four-vector from the mass and the momentum three-vector
    return np.array([np.linalg.norm([m]+p)]+p)

def get_stu(p1,p2,p3):
    def fourproduct(v1,v2):
        return v1[0]*v2[0]-np.dot(v1[1:],v2[1:])
    _s = p1+p2
    _t = p2-p3
    _u = p1-p3
    return fourproduct(_s,_s), fourproduct(_t,_t), fourproduct(_u,_u)

# load c++ library
name = "box1L_ebr"
intlib = IntegralLibrary("{0}/{0}_pylink.so".format(name))
intlib.use_Qmc(transform="korobov3", fitfunction="polysingular", verbosity=1)

# integrate
str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = intlib(real_parameters=[4.0, -2.82842712475, -2.82842712475, 0.1])

# convert complex numbers from c++ to sympy notation
str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')

# convert result to sympy expressions
integral_result = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
integral_result_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))

# print analytic result for given kinematics
p1 = fourvec(0,[1,0,0])
p2 = fourvec(0,[-1,0,0])
p3 = fourvec(0,[0,1,1])
s,t,u = get_stu(p1,p2,p3)
mtsq=0.1

print("When using s, t, u, mtsq = {:.5f}, {:.5f}, {:.5f}, {:.5f}".format(s,t,u,mtsq))
print("The expected result is: {:.15f}".format(f(s,t,u,mtsq)))
print()

# examples how to access individual orders
print('Numerical Result')
for power in [0]:
    valreal, valimg = integral_result.coeff('eps',power).coeff('value').as_real_imag()
    errreal, errimg = integral_result.coeff('eps',power).coeff('error').as_real_imag()
    print("eps^{:<2} {: .15f}{:+.15f}*I +/- {:.15f}{:+.15f}*I".format(power,float(valreal),float(valimg),float(errreal),float(errimg)))
