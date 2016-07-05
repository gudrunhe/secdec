#include <iostream>
#include <stdexcept>
#include <vector>
#include <functional>
#include <cstdint>
#include <cuba.h>

#include "%(name)s.hpp"

// Vegas parameters
#define NCOMP 1
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-16
#define VERBOSE 0
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

struct userdata_t {
    const int epsilon_order;
    %(name)s::real_t const * const real_parameters;
    %(name)s::complex_t const * const complex_parameters;
};

static int cuba_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    userdata_t &data = *(reinterpret_cast<userdata_t *>(userdata));
    ff[0] = 0.;
    for (const auto &sector : %(name)s::sectors)
    {
        if (data.epsilon_order >= sector.get_order_min() && data.epsilon_order <= sector.get_order_max())
            ff[0] += sector[data.epsilon_order].integrand(xx, data.real_parameters, data.complex_parameters);
        // else: does not contribute
    }
    return 0;
}

void integrate(const std::vector<%(name)s::real_t> &real_parameters, const std::vector<%(name)s::complex_t> &complex_parameters)
{

    std::cout << "-- Integrating --" << std::endl;

    if ( real_parameters.size() != %(name)s::number_of_real_parameters )
        throw std::logic_error("Did not set the correct number of real parameters");

    if ( complex_parameters.size() != + %(name)s::number_of_complex_parameters )
        throw std::logic_error("Did not set the correct number of complex parameters");

    int lowest_epsilon_order = %(name)s::sectors.at(0).get_order_min();
    int highest_epsilon_order = %(name)s::sectors.at(0).get_order_max();
    for (const auto &sector : %(name)s::sectors)
    {
        lowest_epsilon_order = std::min(lowest_epsilon_order, sector.get_order_min());
        highest_epsilon_order = std::max(highest_epsilon_order, sector.get_order_max());
    }

    for ( int epsilon_order = lowest_epsilon_order; epsilon_order <= highest_epsilon_order; ++epsilon_order )
    {
        int comp, nregions, neval, fail;
        std::array<cubareal, NCOMP> integral;
        std::array<cubareal, NCOMP> error;
        std::array<cubareal, NCOMP> prob;

        int ndim = 0;
        for (const auto &sector : %(name)s::sectors)
        {
            if (epsilon_order >= sector.get_order_min() && epsilon_order <= sector.get_order_max())
                ndim = std::max(ndim, sector[epsilon_order].number_of_integration_variables);
            // else: does not contribute
        }

        userdata_t userdata{epsilon_order, real_parameters.data(), complex_parameters.data()};

        Vegas(ndim, NCOMP, cuba_integrand, reinterpret_cast<void *>(&userdata), NVEC,
              EPSREL, EPSABS, VERBOSE, SEED,
              MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
              GRIDNO, STATEFILE, SPIN,
              &neval, &fail, integral.data(), error.data(), prob.data()
              );

        std::cout << "epsilon_order: " << epsilon_order << std::endl;
        std::cout << "ndim: " << ndim << std::endl;
        std::cout << "VEGAS RESULT: " << " neval: " << neval << " fail: " << fail << std::endl;

        for( unsigned comp = 0; comp < NCOMP; comp++)
            std::cout << integral.at(comp) << " +- " << error.at(comp) << " p = " << prob.at(comp) << std::endl;

    }
}

int main()
{
    // User Specified Phase-space point
    const std::vector<%(name)s::real_t> real_parameters = { -3., -2. };
    const std::vector<%(name)s::complex_t> complex_parameters = {  };

    integrate(real_parameters, complex_parameters);

    return 0;
}
