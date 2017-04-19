#include <iostream>
#include <complex>
#include <secdecutil/uncertainties.hpp>

int main()
{
    secdecutil::UncorrelatedDeviation<double> r(1.,0.5);
    secdecutil::UncorrelatedDeviation<std::complex<double>> c({2.,3.},{0.6,0.7});

    std::cout << "r: " << r << std::endl;
    std::cout << "c: " << c << std::endl << std::endl;

    std::cout << "r.value:       " << r.value << std::endl;
    std::cout << "r.uncertainty: " << r.uncertainty << std::endl;
    std::cout << "r + c:         " << r + c << std::endl;
    std::cout << "r * c:         " << r * c << std::endl;
    std::cout << "r / 3.0:       " << r / 3. << std::endl;
    // std::cout << "1. / r:     " << 1. / r << std::endl; // ERROR
    // std::cout << "c / r:      " << c / r << std::endl;  // ERROR
}
