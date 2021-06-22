from pySecDec.integral_interface import IntegralLibrary

if __name__ == "__main__":

    # load c++ library
    amp = IntegralLibrary('yyyy1L/yyyy1L_pylink.so')

    # choose Qmc integrator
    amp.use_Qmc()

    # integrate
    _, _, result = amp([1.3,-0.8]) # t, u

    # print result
    print('Numerical Result:' + result)
