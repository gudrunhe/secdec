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
    %(name)s::Mandelstam_t const * const Mandelstam;
    %(name)s::mass_t const * const mass;
};

static int cuba_integrand(const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
    userdata_t &data = *(reinterpret_cast<userdata_t *>(userdata));
    ff[0] = 0.;
    for (const auto &sector : %(name)s::sectors)
    {
        if (data.epsilon_order >= sector.order_min && data.epsilon_order <= sector.order_max)
            ff[0] += sector[data.epsilon_order].integrand(xx, data.Mandelstam, data.mass);
        // else: does not contribute
    }
    return 0;
}

void integrate(const std::vector<%(name)s::Mandelstam_t> &Mandelstam, const std::vector<%(name)s::mass_t> &mass)
{

    std::cout << "-- Integrating --" << std::endl;

    if ( Mandelstam.size() != %(name)s::number_of_Mandelstam_variables )
        throw std::logic_error("Did not set the correct number of Mandelstam variables");

    if ( mass.size() != + %(name)s::number_of_masses )
        throw std::logic_error("Did not set the correct number of masses");

    int lowest_epsilon_order = %(name)s::sectors.at(0).order_min;
    int highest_epsilon_order = %(name)s::sectors.at(0).order_max;
    for (const auto &sector : %(name)s::sectors)
    {
        lowest_epsilon_order = std::min(lowest_epsilon_order, sector.order_min);
        highest_epsilon_order = std::max(highest_epsilon_order, sector.order_max);
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
            if (epsilon_order >= sector.order_min && epsilon_order <= sector.order_max)
                ndim = std::max(ndim, sector[epsilon_order].number_of_Feynman_parameters);
            // else: does not contribute
        }

        userdata_t userdata{epsilon_order, Mandelstam.data(), mass.data()};

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
    const std::vector<%(name)s::Mandelstam_t> Mandelstam = { -2., -1.5 };
    const std::vector<%(name)s::mass_t> mass = { 1. };

    integrate(Mandelstam, mass);

    return 0;
}
