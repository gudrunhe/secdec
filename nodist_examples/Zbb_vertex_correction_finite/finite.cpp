#include <iostream> // std::cout
#include <cmath> // std::log
#include <complex> // std::complex
#include <numeric> // std::accumulate
#include <vector> // std::vector

#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Divonne
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
#include "F1_45_2_alt/F1_45_2_alt.hpp"
#include "F1diminc4_42/F1diminc4_42.hpp"
#include "F1diminc2_37/F1diminc2_37.hpp"
#include "F1diminc2_21/F1diminc2_21.hpp"
#include "F1diminc2_13/F1diminc2_13.hpp"

typedef double real_t;
typedef std::complex<real_t> complex_t;
template<typename T> using nested_series_t = secdecutil::Series<T>;

template<typename integrand_return_t, typename integrand_t, typename ...other_types>
nested_series_t<secdecutil::UncorrelatedDeviation<integrand_return_t>> compute
(
    std::vector<nested_series_t<integrand_t>> integrands,
    nested_series_t<integrand_return_t> prefactor,
    secdecutil::Integrator<integrand_return_t,other_types...> integrator
)
{
    // add integrands of sectors (together flag)
    const nested_series_t<integrand_t> summed_integrands = std::accumulate(++integrands.begin(), integrands.end(), *integrands.begin() );

    // integrate
    return secdecutil::deep_apply(summed_integrands, integrator.integrate) * prefactor;
}


int main()
{
    const std::vector<real_t> real_parameters = {};
    const std::vector<complex_t> complex_parameters = {};

    auto cuhre = secdecutil::cuba::Cuhre<complex_t>();
    cuhre.flags = 2; // verbose output
    cuhre.epsrel = 1e-8;
    cuhre.epsabs = 1e-8;
    cuhre.maxeval = 1e6;

    auto divonne = secdecutil::cuba::Divonne<complex_t>();
    divonne.flags = 2; // verbose output
    divonne.epsrel = 1e-8;
    divonne.epsabs = 1e-8;
    divonne.maxeval = 1e6;
    divonne.border = 1e-8;


    // F1diminc2_37::nested_series_t<secdecutil::UncorrelatedDeviation<F1diminc2_37::integrand_return_t>> f1diminc2x37xRes = f1diminc2x37(); // factorizable - requires split
    // Insert analytic result for factorizable F1diminc2_37
    const secdecutil::Series<complex_t> f1diminc2x37xRes =
    {
        0,4,
        {
            {0.166666666666666666666666666667,0.},
            {0.333333333333333333333333333334,0.523598775598298873077107230546},
            {-0.155800366757446551569540916655,1.047197551196597746154214461093},
            {-0.84584824603026834109452101620,1.23310963905153382076742017032},
            {-1.55640608784664440824954157034,0.78783121758614367388885109920}
        },
        true,"eps"
    };

    // F1_45_2::nested_series_t<secdecutil::UncorrelatedDeviation<F1_45_2::integrand_retucd rn_t>> f1x45x2xRes = f1x45x2(); // requires split

    #define COMPUTE(RESULT, NAME, DEFORMATION_PARAMETERS_MAXIMUM, INTEGRATOR) \
    \
    NAME::nested_series_t<secdecutil::UncorrelatedDeviation<NAME::integrand_return_t>> RESULT = \
    compute \
    ( \
        NAME::make_integrands( \
                                        real_parameters, \
                                        complex_parameters, \
                                        10000, /*number of presamples*/ \
                                        DEFORMATION_PARAMETERS_MAXIMUM \
                                    ), \
        NAME::prefactor(real_parameters,complex_parameters), \
        INTEGRATOR \
    );

    COMPUTE(f1x45x2xAltxRes, F1_45_2_alt, 1., divonne);
    COMPUTE(f1x45xRes, F1_45, 1., divonne);
    COMPUTE(f1diminc2x63xRes, F1diminc2_63, 0.1, divonne);

    COMPUTE(f1diminc2x62xRes, F1diminc2_62, 0.1, divonne);
    COMPUTE(f1diminc2x61xRes, F1diminc2_61, 0.1, divonne);
    COMPUTE(f1x47xRes, F1_47, 0.1, divonne);

    COMPUTE(f1diminc4x51xRes, F1diminc4_51, 0.1, divonne);
    COMPUTE(f1diminc2x46xRes, F1diminc2_46, 1., divonne);
    COMPUTE(f1diminc4x42xRes, F1diminc4_42, 1., cuhre);

    COMPUTE(f1diminc2x21xRes, F1diminc2_21, 1., cuhre);
    COMPUTE(f1diminc2x13xRes, F1diminc2_13, 1., cuhre);


    std::cout << "f1diminc2x63xRes " << f1diminc2x63xRes << std::endl;
    std::cout << "f1diminc2x62xRes " << f1diminc2x62xRes << std::endl;
    std::cout << "f1diminc2x61xRes " << f1diminc2x61xRes << std::endl;
    std::cout << "f1x47xRes " << f1x47xRes << std::endl;
    std::cout << "f1diminc4x51xRes " << f1diminc4x51xRes << std::endl;
    std::cout << "f1diminc2x46xRes " << f1diminc2x46xRes << std::endl;
    std::cout << "f1x45xRes " << f1x45xRes << std::endl;
    //    std::cout << "f1x45x2xRes " << f1x45x2xRes << std::endl;
    std::cout << "f1x45x2xAltxRes " << f1x45x2xAltxRes << std::endl;
    std::cout << "f1diminc2x37xRes " << f1diminc2x37xRes << std::endl;
    std::cout << "f1diminc4x42xRes " << f1diminc4x42xRes << std::endl;
    std::cout << "f1diminc2x21xRes " << f1diminc2x21xRes << std::endl;
    std::cout << "f1diminc2x13xRes " << f1diminc2x13xRes << std::endl;

    const secdecutil::Series<real_t> eps = {1,1,{1},false,"eps"};

    // Finite basis 1
    //    auto result =
    //    (
    //     + 1./(eps*eps*eps*eps)*(-2*f1diminc2x13xRes + f1diminc2x21xRes + 4*f1diminc4x42xRes - 2*f1diminc4x51xRes)/2.
    //     + 1./(eps*eps*eps)*(2*f1diminc2x13xRes - 15*f1diminc2x21xRes + 18*f1diminc2x37xRes - 68*f1diminc4x42xRes + 38*f1diminc4x51xRes)/4.
    //     + 1./(eps*eps)*(-14*f1diminc2x13xRes + 43*f1diminc2x21xRes - 27*f1diminc2x37xRes - 4*f1diminc2x46xRes + 92*f1diminc4x42xRes - 120*f1diminc4x51xRes + 8*f1x45x2xRes)/4.
    //     + 1./(eps)*(-3*f1diminc2x13xRes - 25*f1diminc2x21xRes - 57*f1diminc2x37xRes + 8*f1diminc2x46xRes + 2*f1diminc2x61xRes + f1diminc2x62xRes + 2*f1diminc2x63xRes + 44*f1diminc4x42xRes + 108*f1diminc4x51xRes + 8*f1x45x2xRes + 9*f1x45xRes)/4.
    //     + (5*f1diminc2x13xRes - 41*f1diminc2x21xRes - 141*f1diminc2x37xRes + 14*f1diminc2x46xRes - 6*f1diminc2x61xRes - 3*f1diminc2x62xRes -
    //        10*f1diminc2x63xRes + 12*f1diminc4x42xRes + 8*f1x45x2xRes + 6*f1x45xRes + f1x47xRes)/4.
    //     );

    // Finite basis 2 - avoids f1x45x2 which has squared singular F poly in denominator
    auto result =
    + 1./(eps*eps*eps*eps)*(f1diminc2x13xRes + f1diminc2x21xRes/2. - 3.*f1diminc2x37xRes + 2.*f1diminc4x42xRes - f1diminc4x51xRes)
    + 1./(eps*eps*eps)*((-18.*f1diminc2x13xRes - 15.*f1diminc2x21xRes + 2.*(15.*f1diminc2x37xRes - 34.*f1diminc4x42xRes + 19.*f1diminc4x51xRes + f1x45x2xAltxRes + f1x45xRes))/4.)
    + 1./(eps*eps)*((2.*f1diminc2x13xRes + 43.*f1diminc2x21xRes - 15.*f1diminc2x37xRes - 4.*f1diminc2x46xRes + 92.*f1diminc4x42xRes - 120.*f1diminc4x51xRes -
                     8.*f1x45x2xAltxRes - 6.*f1x45xRes)/4.)
    + 1./(eps)*((-15.*f1diminc2x13xRes - 25.*f1diminc2x21xRes - 21.*f1diminc2x37xRes + 8.*f1diminc2x46xRes + 2.*f1diminc2x61xRes + f1diminc2x62xRes +
                 2.*f1diminc2x63xRes + 44.*f1diminc4x42xRes + 108.*f1diminc4x51xRes + 4.*f1x45x2xAltxRes + f1x45xRes)/4.)
    + ((-31.*f1diminc2x13xRes - 41.*f1diminc2x21xRes - 33.*f1diminc2x37xRes + 14.*f1diminc2x46xRes - 6.*f1diminc2x61xRes - 3.*f1diminc2x62xRes -
        10.*f1diminc2x63xRes + 12.*f1diminc4x42xRes - 6.*f1x45xRes + f1x47xRes)/4.)
    ;
    
    std::cout << result << std::endl;
    
    return 0;
}
