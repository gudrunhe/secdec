#include <iostream>
#include <complex>
#include <secdecutil/series.hpp>
#include <secdecutil/deep_apply.hpp>

int main()
{
    std::function<std::complex<double>(std::complex<double>)> conjugate =
    [] (std::complex<double> element)
    {
        return std::conj(element);
    };

    secdecutil::Series<std::complex<double>> u(-1,0,{{1,2},{3,4}},false,"eps");
    secdecutil::Series<secdecutil::Series<std::complex<double>>> m(1,1,{{1,1,{{1,2}},false,"alpha"},},false,"eps");

    std::cout << "u: " << u << std::endl;
    std::cout << "m: " << m << std::endl << std::endl;

    std::cout << "conjugated u:   " << secdecutil::deep_apply(u, conjugate) << std::endl;
    std::cout << "conjugated m: " << secdecutil::deep_apply(m, conjugate) << std::endl;
}
