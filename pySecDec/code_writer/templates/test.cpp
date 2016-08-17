#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <numeric> // std::accumulate
#include <functional> // std::bind
#include <type_traits> // std::remove_const
#include <typeinfo>

#include <secdecutil/series.hpp> // Series
#include <secdecutil/uncertainties.hpp> // GaussianUncertainty
#include <secdecutil/sector_container.hpp> // SectorContainer to IntegrandContainer
#include <secdecutil/integrand_container.hpp> // IntegrandContainer
#include <secdecutil/deep_apply.hpp> // deep_apply

#include "%(name)s.hpp"

void print_integral_info()
{
    std::cout << "-- print_integral_info --" << std::endl;
    std::cout << "%(name)s::number_of_sectors " << %(name)s::number_of_sectors << std::endl;
    std::cout << "%(name)s::number_of_regulators " << %(name)s::number_of_regulators << std::endl;
    std::cout << "%(name)s::number_of_real_parameters " << %(name)s::number_of_real_parameters << std::endl;
    std::cout << "%(name)s::number_of_complex_parameters " << %(name)s::number_of_complex_parameters << std::endl;

    std::cout << "%(name)s::lowest_orders";
    for ( const auto& lowest_order : %(name)s::lowest_orders )
        std::cout << " " << lowest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::highest_orders";
    for ( const auto& highest_order : %(name)s::highest_orders )
        std::cout << " " << highest_order;
    std::cout << std::endl;

    std::cout << "%(name)s::requested_orders";
    for ( const auto& requested_order : %(name)s::requested_orders )
        std::cout << " " << requested_order;
    std::cout << std::endl;
}

namespace cuba //TODO: move to SecDecUtil
{
    #include <cuba.h>

    template<typename T>
    struct CubaVegas
    {
        static int cuba_integrand_prototype(const int *ndim, const T integration_variables[], const int *ncomp, T result[], void *userdata)
        {
            auto& integrand_container = *( reinterpret_cast<secdecutil::IntegrandContainer<T, T const * const> *>(userdata) );
            result[0] = integrand_container.integrand(integration_variables);
            return 0;
        };
        constexpr static int ncomp = 1;
        constexpr static int nvec = 1;
        T epsrel;
        T epsabs;
        int flags;
        int seed;
        int mineval;
        int maxeval;
        int nstart;
        int nincrease;
        int nbatch;
        constexpr static int gridno = 0;
        constexpr static char * statefile = nullptr;
        constexpr static void* spin = nullptr;

        CubaVegas(
            T epsrel = 1e-2,
            T epsabs = 1e-16,
            int flags = 0,
            int seed = 0,
            int mineval = 0,
            int maxeval = 50000,
            int nstart = 1000,
            int nincrease = 500,
            int nbatch = 1000
        ) :
            epsrel(epsrel),epsabs(epsabs),
            flags(flags),seed(seed),mineval(mineval),maxeval(maxeval),
            nstart(nstart),nincrease(nincrease),nbatch(nbatch)
        {};

        std::function<secdecutil::GaussianUncertainty<T>(secdecutil::IntegrandContainer<T, T const * const>)>
        integrate =
        [ this ] (secdecutil::IntegrandContainer<T, T const * const> integrand_container)
        {
            std::cout << "-- Integrating --" << std::endl;

            // Cuba output values
            std::array<T, ncomp> integral;
            std::array<T, ncomp> error;
            std::array<T, ncomp> prob;
            int comp, nregions, neval, fail;

            // Cuba call
            Vegas(
                  integrand_container.number_of_integration_variables,
                  ncomp,
                  cuba_integrand_prototype,
                  reinterpret_cast<void*>(&integrand_container), // userdata
                  nvec,
                  epsrel,
                  epsabs,
                  flags,
                  seed,
                  mineval,
                  maxeval,
                  nstart,
                  nincrease,
                  nbatch,
                  gridno,
                  statefile,
                  spin,
                  &neval,
                  &fail,
                  integral.data(),
                  error.data(),
                  prob.data()
                  );

            std::cout << "VEGAS RESULT: " << " neval: " << neval << " fail: " << fail << std::endl;

            for( unsigned comp = 0; comp < ncomp; comp++)
                std::cout << integral.at(comp) << " +- " << error.at(comp) << " p = " << prob.at(comp) << std::endl;

            return secdecutil::GaussianUncertainty<T>(integral.at(0),error.at(0));
        };
    };
    template<typename T>
    struct CubaVegas<std::complex<T>>
    {
        static int cuba_integrand_prototype(const int *ndim, const T integration_variables[], const int *ncomp, T result[], void *userdata)
        {
            auto& integrand_container = *( reinterpret_cast<secdecutil::IntegrandContainer<std::complex<T>, T const * const> *>(userdata) );
            std::complex<T> evaluated_integrand = integrand_container.integrand(integration_variables);
            result[0] = evaluated_integrand.real();
            result[1] = evaluated_integrand.imag();
            return 0;
        };
        constexpr static int ncomp = 2;
        constexpr static int nvec = 1;
        T epsrel;
        T epsabs;
        int flags;
        int seed;
        int mineval;
        int maxeval;
        int nstart;
        int nincrease;
        int nbatch;
        constexpr static int gridno = 0;
        constexpr static char * statefile = nullptr;
        constexpr static void* spin = nullptr;

        std::function<secdecutil::GaussianUncertainty<std::complex<T>>(secdecutil::IntegrandContainer<std::complex<T>, T const * const>)>
        integrate =
        [ this ] (secdecutil::IntegrandContainer<std::complex<T>, T const * const> integrand_container)
        {
            std::cout << "-- Integrating --" << std::endl;

            // Cuba output values
            std::array<T, ncomp> integral;
            std::array<T, ncomp> error;
            std::array<T, ncomp> prob;
            int comp, nregions, neval, fail;

            // Cuba call
            Vegas(
                  integrand_container.number_of_integration_variables,
                  ncomp,
                  cuba_integrand_prototype,
                  reinterpret_cast<void*>(&integrand_container), // userdata
                  nvec,
                  epsrel,
                  epsabs,
                  flags,
                  seed,
                  mineval,
                  maxeval,
                  nstart,
                  nincrease,
                  nbatch,
                  gridno,
                  statefile,
                  spin,
                  &neval,
                  &fail,
                  integral.data(),
                  error.data(),
                  prob.data()
                  );

            std::cout << "VEGAS RESULT: " << " neval: " << neval << " fail: " << fail << std::endl;

            for( unsigned comp = 0; comp < ncomp; comp++)
                std::cout << integral.at(comp) << " +- " << error.at(comp) << " p = " << prob.at(comp) << std::endl;

            return secdecutil::GaussianUncertainty<std::complex<T>>({integral.at(0),integral.at(1)},{error.at(0),error.at(1)});
        };

        CubaVegas(
            T epsrel = 1e-2,
            T epsabs = 1e-16,
            int flags = 0,
            int seed = 0,
            int mineval = 0,
            int maxeval = 50000,
            int nstart = 1000,
            int nincrease = 500,
            int nbatch = 1000
        ) :
            epsrel(epsrel),epsabs(epsabs),
            flags(flags),seed(seed),mineval(mineval),maxeval(maxeval),
            nstart(nstart),nincrease(nincrease),nbatch(nbatch)
        {};
    };
};

int main()
{
    // TODO - write method to parse arguments and check validity
    // User Specified Phase-space point
    const std::vector<%(name)s::real_t> real_parameters = { -3., -2. };
    const std::vector<%(name)s::complex_t> complex_parameters = {  };
    if ( real_parameters.size() != %(name)s::number_of_real_parameters )
        throw std::logic_error("Did not set the correct number of real parameters");
    if ( complex_parameters.size() != + %(name)s::number_of_complex_parameters )
        throw std::logic_error("Did not set the correct number of complex parameters");

    print_integral_info();

    const auto sector_integrands = %(name)s::make_integrands(real_parameters, complex_parameters);

    // const auto integrand_container_sum = sector_integrands.at(0) + sector_integrands.at(1); // Example how to add integrand containers

    // Add integrands of sectors (together flag)
    const auto all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

    // Integrate
    auto integrator = cuba::CubaVegas<%(name)s::integrand_return_t>();
    auto result_all = secdecutil::deep_apply( all_sectors,  integrator.integrate );

    std::cout << "-- All -- " << std::endl;
    std::cout << result_all << std::endl;

    return 0;
}
