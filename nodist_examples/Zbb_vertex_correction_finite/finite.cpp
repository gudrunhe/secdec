#include <iostream> // std::cout
#include <cmath> // std::log
#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas
#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include "F1diminc2_63/F1diminc2_63.hpp"
#include "F1diminc2_62/F1diminc2_62.hpp"
#include "F1diminc2_61/F1diminc2_61.hpp"
#include "F1_47/F1_47.hpp"
#include "F1diminc4_51/F1diminc4_51.hpp"
#include "F1diminc2_46/F1diminc2_46.hpp"
#include "F1_45/F1_45.hpp"
#include "F1_45_2/F1_45_2.hpp"
#include "F1diminc4_42/F1diminc4_42.hpp"
#include "F1diminc2_37/F1diminc2_37.hpp"
#include "F1diminc2_21/F1diminc2_21.hpp"
#include "F1diminc2_13/F1diminc2_13.hpp"

/*
 * pySecDec Master Integrals
 */
F1diminc2_63::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_63::integrand_return_t>> f1diminc2x63()
{
    using namespace F1diminc2_63;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_63::integrand_t>> integrands = F1diminc2_63::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_63::nested_series_t<F1diminc2_63::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_63::prefactor(real_parameters, complex_parameters);
}

F1diminc2_62::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_62::integrand_return_t>> f1diminc2x62()
{
    using namespace F1diminc2_62;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_62::integrand_t>> integrands = F1diminc2_62::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_62::nested_series_t<F1diminc2_62::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_62::prefactor(real_parameters, complex_parameters);
}

F1diminc2_61::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_61::integrand_return_t>> f1diminc2x61()
{
    using namespace F1diminc2_61;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_61::integrand_t>> integrands = F1diminc2_61::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_61::nested_series_t<F1diminc2_61::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_61::prefactor(real_parameters, complex_parameters);
}

F1_47::nested_series_t<secdecutil::UncorrelatedDeviation<F1_47::integrand_return_t>> f1x47()
{
    using namespace F1_47;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1_47::integrand_t>> integrands = F1_47::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                               // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                               // optionally set here as additional arguments.
                                                                                               );

    // add integrands of sectors (together flag)
    const F1_47::nested_series_t<F1_47::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1_47::prefactor(real_parameters, complex_parameters);
}

F1diminc4_51::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc4_51::integrand_return_t>> f1diminc4x51()
{
    using namespace F1diminc4_51;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc4_51::integrand_t>> integrands = F1diminc4_51::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc4_51::nested_series_t<F1diminc4_51::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc4_51::prefactor(real_parameters, complex_parameters);
}

F1diminc2_46::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_46::integrand_return_t>> f1diminc2x46()
{
    using namespace F1diminc2_46;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_46::integrand_t>> integrands = F1diminc2_46::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_46::nested_series_t<F1diminc2_46::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_46::prefactor(real_parameters, complex_parameters);
}

F1_45::nested_series_t<secdecutil::UncorrelatedDeviation<F1_45::integrand_return_t>> f1x45()
{
    using namespace F1_45;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1_45::integrand_t>> integrands = F1_45::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                               // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                               // optionally set here as additional arguments.
                                                                                               );

    // add integrands of sectors (together flag)
    const F1_45::nested_series_t<F1_45::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1_45::prefactor(real_parameters, complex_parameters);
}

F1_45_2::nested_series_t<secdecutil::UncorrelatedDeviation<F1_45_2::integrand_return_t>> f1x45x2()
{
    using namespace F1_45_2;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1_45_2::integrand_t>> integrands = F1_45_2::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                   // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                   // optionally set here as additional arguments.
                                                                                                   );

    // add integrands of sectors (together flag)
    const F1_45_2::nested_series_t<F1_45_2::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1_45_2::prefactor(real_parameters, complex_parameters);
}

F1diminc4_42::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc4_42::integrand_return_t>> f1diminc4x42()
{
    using namespace F1diminc4_42;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc4_42::integrand_t>> integrands = F1diminc4_42::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc4_42::nested_series_t<F1diminc4_42::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc4_42::prefactor(real_parameters, complex_parameters);
}

F1diminc2_37::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_37::integrand_return_t>> f1diminc2x37()
{
    using namespace F1diminc2_37;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 1e2;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_37::integrand_t>> integrands = F1diminc2_37::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_37::nested_series_t<F1diminc2_37::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_37::prefactor(real_parameters, complex_parameters);
}

F1diminc2_21::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_21::integrand_return_t>> f1diminc2x21()
{
    using namespace F1diminc2_21;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_21::integrand_t>> integrands = F1diminc2_21::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_21::nested_series_t<F1diminc2_21::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_21::prefactor(real_parameters, complex_parameters);
}

F1diminc2_13::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_13::integrand_return_t>> f1diminc2x13()
{
    using namespace F1diminc2_13;

    const std::vector<real_t> real_parameters{};
    const std::vector<complex_t> complex_parameters{};
    const unsigned number_of_samples = 100000;
    const double deformation_parameters_maximum = 0.1;

    // optimize contour
    const std::vector<nested_series_t<F1diminc2_13::integrand_t>> integrands = F1diminc2_13::make_integrands(real_parameters, complex_parameters, number_of_samples, deformation_parameters_maximum
                                                                                                             // The number of samples for the contour optimization, the minimal and maximal deformation parameters, and the decrease factor can be
                                                                                                             // optionally set here as additional arguments.
                                                                                                             );

    // add integrands of sectors (together flag)
    const F1diminc2_13::nested_series_t<F1diminc2_13::integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // define the integrator
    auto integrator = secdecutil::cuba::Vegas<std::complex<double>>();
    integrator.flags = 2; // verbose output
    integrator.epsrel = 1e-5;
    integrator.epsabs = 1e-7;
    integrator.maxeval = 1e6;

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * F1diminc2_13::prefactor(real_parameters, complex_parameters);
}

int main()
{

    F1diminc2_37::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_37::integrand_return_t>> f1diminc2x37xRes = f1diminc2x37(); // factorizable - requires split
    F1_45_2::nested_series_t<secdecutil::UncorrelatedDeviation<F1_45_2::integrand_return_t>> f1x45x2xRes = f1x45x2(); // requires split
    F1_45::nested_series_t<secdecutil::UncorrelatedDeviation<F1_45::integrand_return_t>> f1x45xRes = f1x45();
    F1diminc2_63::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_63::integrand_return_t>> f1diminc2x63xRes = f1diminc2x63();
    F1diminc2_62::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_62::integrand_return_t>> f1diminc2x62xRes = f1diminc2x62();
    F1diminc2_61::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_61::integrand_return_t>> f1diminc2x61xRes = f1diminc2x61();
    F1_47::nested_series_t<secdecutil::UncorrelatedDeviation<F1_47::integrand_return_t>> f1x47xRes = f1x47();
    F1diminc4_51::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc4_51::integrand_return_t>> f1diminc4x51xRes = f1diminc4x51();
    F1diminc2_46::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_46::integrand_return_t>> f1diminc2x46xRes = f1diminc2x46();
    F1diminc4_42::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc4_42::integrand_return_t>> f1diminc4x42xRes = f1diminc4x42();
    F1diminc2_21::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_21::integrand_return_t>> f1diminc2x21xRes = f1diminc2x21();
    F1diminc2_13::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_13::integrand_return_t>> f1diminc2x13xRes = f1diminc2x13();

    std::cout << "f1diminc2x63xRes " << f1diminc2x63xRes << std::endl;
    std::cout << "f1diminc2x62xRes " << f1diminc2x62xRes << std::endl;
    std::cout << "f1diminc2x61xRes " << f1diminc2x61xRes << std::endl;
    std::cout << "f1x47xRes " << f1x47xRes << std::endl;
    std::cout << "f1diminc4x51xRes " << f1diminc4x51xRes << std::endl;
    std::cout << "f1diminc2x46xRes " << f1diminc2x46xRes << std::endl;
    std::cout << "f1x45xRes " << f1x45xRes << std::endl;
    std::cout << "f1x45x2xRes " << f1x45x2xRes << std::endl;
    std::cout << "f1diminc2x37xRes " << f1diminc2x37xRes << std::endl;
    std::cout << "f1diminc4x42xRes " << f1diminc4x42xRes << std::endl;
    std::cout << "f1diminc2x21xRes " << f1diminc2x21xRes << std::endl;
    std::cout << "f1diminc2x13xRes " << f1diminc2x13xRes << std::endl;

    const secdecutil::Series<double> eps = {1,1,{1},false,"eps"};

    auto result =
    (
     + 1./(eps*eps*eps*eps)*(-2*f1diminc2x13xRes + f1diminc2x21xRes + 4*f1diminc4x42xRes - 2*f1diminc4x51xRes)/2.
     + 1./(eps*eps*eps)*(2*f1diminc2x13xRes - 15*f1diminc2x21xRes + 18*f1diminc2x37xRes - 68*f1diminc4x42xRes + 38*f1diminc4x51xRes)/4.
     + 1./(eps*eps)*(-14*f1diminc2x13xRes + 43*f1diminc2x21xRes - 27*f1diminc2x37xRes - 4*f1diminc2x46xRes + 92*f1diminc4x42xRes - 120*f1diminc4x51xRes + 8*f1x45x2xRes)/4.
     + 1./(eps)*(-3*f1diminc2x13xRes - 25*f1diminc2x21xRes - 57*f1diminc2x37xRes + 8*f1diminc2x46xRes + 2*f1diminc2x61xRes + f1diminc2x62xRes + 2*f1diminc2x63xRes + 44*f1diminc4x42xRes + 108*f1diminc4x51xRes + 8*f1x45x2xRes + 9*f1x45xRes)/4.
     + (5*f1diminc2x13xRes - 41*f1diminc2x21xRes - 141*f1diminc2x37xRes + 14*f1diminc2x46xRes - 6*f1diminc2x61xRes - 3*f1diminc2x62xRes -
        10*f1diminc2x63xRes + 12*f1diminc4x42xRes + 8*f1x45x2xRes + 6*f1x45xRes + f1x47xRes)/4.
     );

    std::cout << result << std::endl;

    return 0;
}
