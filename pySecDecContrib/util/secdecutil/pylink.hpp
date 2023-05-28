#ifndef SecDecUtil_pylink_hpp_included
#define SecDecUtil_pylink_hpp_included

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <secdecutil/deep_apply.hpp> // deep_apply
#include <secdecutil/integrators/cquad.hpp> // CQuad
#include <secdecutil/integrators/qmc.hpp> // Qmc
#include <secdecutil/integrators/cuba.hpp> // Vegas, Suave, Divonne, Cuhre

#if integral_need_complex 
    #define SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX(NOARGS) do {integrator->together = real_complex_together;} while (false)
#else
    #define SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX(NOARGS)
#endif

using namespace INTEGRAL_NAME;
using secdec_integrand_t = INTEGRAL_NAME::integrand_t; // avoid name conflict with cuba
#ifdef SECDEC_WITH_CUDA
    using INTEGRAL_NAME::cuda_integrand_t;
    using INTEGRAL_NAME::cuda_together_integrand_t;
#endif

extern "C"
{
    /*
     * string (de)allocation
     */
    const char * string2charptr(std::string * str)
    {
        return str->c_str();
    }
    std::string * allocate_string()
    {
        return new std::string;
    }
    void free_string(std::string * strptr)
    {
        delete strptr;
    }


    /*
     * integrator (de)allocation
     */
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_MultiIntegrator(
                                 secdecutil::Integrator<integrand_return_t,real_t>* low_dim_integrator,
                                 secdecutil::Integrator<integrand_return_t,real_t>* high_dim_integrator,
                                 int critical_dim
                            )
    {
        auto integrator = new secdecutil::MultiIntegrator<integrand_return_t,real_t>
            (*low_dim_integrator,*high_dim_integrator,critical_dim);
        return integrator;
    }

    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_gsl_cquad(
                            double epsrel,
                            double epsabs,
                            unsigned int n,
                            bool verbose,
                            double zero_border
                       )
    {
        auto integrator = new secdecutil::gsl::CQuad<integrand_return_t>
            (epsrel,epsabs,n,verbose,zero_border);
        return integrator;
    }

    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Vegas(
                            double epsrel,
                            double epsabs,
                            int flags,
                            int seed,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            long long int nstart,
                            long long int nincrease,
                            long long int nbatch,
                            bool real_complex_together
                       )
    {
        auto integrator = new secdecutil::cuba::Vegas<integrand_return_t>
            (epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,nstart,nincrease,nbatch);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Suave(
                            double epsrel,
                            double epsabs,
                            int flags,
                            int seed,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            long long int nnew,
                            long long int nmin,
                            double flatness,
                            bool real_complex_together
                       )
    {
        auto integrator = new secdecutil::cuba::Suave<integrand_return_t>
            (epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,nnew,nmin,flatness);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Divonne(
                            double epsrel,
                            double epsabs,
                            int flags,
                            int seed,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            int key1,
                            int key2,
                            int key3,
                            int maxpass,
                            double border,
                            double maxchisq,
                            double mindeviation,
                            bool real_complex_together
                         )
    {
        auto integrator = new secdecutil::cuba::Divonne<integrand_return_t>
            (epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,key1,key2,key3,maxpass,
             border,maxchisq,mindeviation);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    secdecutil::Integrator<integrand_return_t,real_t> *
    allocate_cuba_Cuhre(
                            double epsrel,
                            double epsabs,
                            int flags,
                            long long int mineval,
                            long long int maxeval,
                            double zero_border,
                            int key,
                            bool real_complex_together
                       )
    {
        auto integrator = new secdecutil::cuba::Cuhre<integrand_return_t>
            (epsrel,epsabs,flags,mineval,maxeval,zero_border,key);
        SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX();
        return integrator;
    }
    // qmc allocate function prototypes (implementation is done in pylink.cpp)
    #define COMMON_ALLOCATE_QMC_ARGS \
    double epsrel, \
    double epsabs, \
    unsigned long long int maxeval, \
    int errormode, \
    unsigned long long int evaluateminn, \
    unsigned long long int minn, \
    unsigned long long int minm, \
    unsigned long long int maxnperpackage, \
    unsigned long long int maxmperpackage, \
    long long int cputhreads, \
    unsigned long long int cudablocks, \
    unsigned long long int cudathreadsperblock, \
    unsigned long long int verbosity, \
    long long int seed, \
    int transform_id, \
    int fitfunction_id, \
    int generatingvectors_id, \
    unsigned long long int lattice_candidates, \
    bool standard_lattices, \
    bool keep_lattices
    #ifdef SECDEC_WITH_CUDA
        secdecutil::Integrator<integrand_return_t,real_t,cuda_together_integrand_t> *
        allocate_cuda_integrators_Qmc_together(
                                                   COMMON_ALLOCATE_QMC_ARGS,
                                                   unsigned long long int number_of_devices,
                                                   int devices[]
                                              );
        secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t> *
        allocate_cuda_integrators_Qmc_separate(
                                                   COMMON_ALLOCATE_QMC_ARGS,
                                                   unsigned long long int number_of_devices,
                                                   int devices[]
                                              );
    #else
        secdecutil::Integrator<integrand_return_t,real_t> * allocate_integrators_Qmc(COMMON_ALLOCATE_QMC_ARGS);
    #endif
    #undef COMMON_ALLOCATE_QMC_ARGS
    void free_integrator (secdecutil::Integrator<integrand_return_t,real_t> * integrator)
    {
        delete integrator;
    }
    #ifdef SECDEC_WITH_CUDA
        void free_cuda_together_integrator (secdecutil::Integrator<integrand_return_t,real_t,cuda_together_integrand_t> * integrator)
        {
            delete integrator;
        }
        void free_cuda_separate_integrator (secdecutil::Integrator<integrand_return_t,real_t,cuda_integrand_t> * integrator)
        {
            delete integrator;
        }
    #endif

    #undef SET_INTEGRATOR_TOGETHER_OPTION_IF_COMPLEX

}

#endif
