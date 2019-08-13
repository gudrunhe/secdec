from .amplitude_interface import *
import unittest
from nose.plugins.attrib import attr
import re
import os.path, sys, shutil

python_major_version = sys.version[0]

class TestConfiguration(unittest.TestCase):

    #@attr('active')
    def test_configure(self):
        integrator1 = CubaIntegrator("Vegas",epsrel="1e-03", epsabs=1e-5)
        integrator2 = Qmc(transform="korobov4x90")
        contourdef1 = ContourDeformation(number_of_presamples=33000)
        contourdef2 = ContourDeformation()

        testfolder =  'tmpdir_test_configure_python' + python_major_version
        os.mkdir(testfolder)

        try:
            with open(os.path.join(testfolder,"integral_names.txt"),"w") as f:
                f.write("integral1\n\nsubintegral1\nsubintegral2")

            config_filename = os.path.join(testfolder,"config_integral1.hpp")

            conf = Configuration(testfolder)
            conf.set_options(integrator1, contourdef1)
            conf.set_options(integrator2, contourdef2, sub_integral_name="subintegral2")

            conf.configure()

            target_configure_file = """#ifndef integral1_config_integral1_hpp_included
#define integral1_config_integral1_hpp_included

/*
 * Set contour-deformation and integrator properties. These can be set per
 * integral by defining an Options struct in the appropriate namespace
 * (see below). Integranls for which no special Options struct is defined
 * use the Options struct in the parent namespace of all the integrands.
 */

namespace integral1
{
    // This struct contains the options for the integrals of the enclosing namespace.
    struct Options
    {
        // integrator options
        // whether the integrator to be used supports CUDA (out of the builtin integrators, only the Qmc does)
        constexpr static bool cuda_compliant_integrator = false;

        // integral settings
        template<typename integrand_t, unsigned long long int maximal_number_of_integration_variables>
        struct Integral
        {
            // the integrator (Cuhre, Divonne, Qmc, Suave, Vegas)
            using integrator_t = secdecutil::cuba::Vegas<integrand_return_t>;

            // instantiation and options of the integrator
            static std::shared_ptr<integrator_t> configure_integrator()
            {
                std::shared_ptr<integrator_t> integrator = std::make_shared<integrator_t>();
                integrator->epsrel = 1e-03;
                integrator->epsabs = 1e-05;
                return integrator;
            }

            // the kind of integral (CubaIntegral, QmcIntegral, <CustomIntegral>)
            using integral_t = secdecutil::amplitude::CubaIntegral<integrand_return_t,real_t,integrator_t,integrand_t>;
        };

        // contour deformation parameters
        struct ContourDeformation
        {
            constexpr static unsigned number_of_presamples = 33000;
            constexpr static real_t deformation_parameters_maximum = 1.0;
            constexpr static real_t deformation_parameters_minimum = 1e-05;
            constexpr static real_t deformation_parameters_decrease_factor = 0.9;
        };
    };

    namespace subintegral2
    {
        // This struct contains the options for the integrals of the enclosing namespace.
        struct Options
        {
            // integrator options
            // whether the integrator to be used supports CUDA (out of the builtin integrators, only the Qmc does)
            constexpr static bool cuda_compliant_integrator = true;

            // integral settings
            template<typename integrand_t, unsigned long long int maximal_number_of_integration_variables>
            struct Integral
            {
                // the integrator (Cuhre, Divonne, Qmc, Suave, Vegas)
                using integrator_t = secdecutil::integrators::Qmc<
                                                                        integrand_return_t,
                                                                        maximal_number_of_integration_variables,
                                                                        ::integrators::transforms::Korobov<4,90>::type,
                                                                        integrand_t
                                                                   >;

                // instantiation and options of the integrator
                static std::shared_ptr<integrator_t> configure_integrator()
                {
                    std::shared_ptr<integrator_t> integrator = std::make_shared<integrator_t>();
                    // no integrator options set
                    return integrator;
                }

                // the kind of integral (CubaIntegral, QmcIntegral, <CustomIntegral>)
                using integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,integrator_t,integrand_t>;
            };

            // contour deformation parameters
            struct ContourDeformation
            {
                constexpr static unsigned number_of_presamples = 100000;
                constexpr static real_t deformation_parameters_maximum = 1.0;
                constexpr static real_t deformation_parameters_minimum = 1e-05;
                constexpr static real_t deformation_parameters_decrease_factor = 0.9;
            };
        };
    }
}

#endif
"""

            with open(config_filename) as file:
                contents = file.read()

            # sort the integrator options alphabetically, needed for python version < 3.6, because there dicts are unordered
            def sort_options(contents):
                return re.sub(r"(?<=std::shared_ptr<integrator_t> integrator = std::make_shared<integrator_t>\(\);\n)(.*?)(?=\s*return integrator;)",
                    lambda x: "\n".join(sorted(x.group(0).splitlines())), contents, flags=re.DOTALL)

            contents = sort_options(contents)
            target_configure_file = sort_options(target_configure_file)

            self.assertEqual(contents, target_configure_file)


        finally:
            shutil.rmtree(testfolder)

class TestQmc(unittest.TestCase):

    #@attr('active')
    def test_transform(self):

        # no default transform
        self.assertRaisesRegexp(TypeError, '(M|m)issing.*transform', Qmc)

        # check some transforms
        transforms_to_test = ["none", "baker", "korobov1x1", "Korobov2x5", "kOrObOv4", 'Sidi8']
        target_replacements = [
                                   '::integrators::transforms::None::type',
                                   '::integrators::transforms::Baker::type',
                                   '::integrators::transforms::Korobov<1,1>::type',
                                   '::integrators::transforms::Korobov<2,5>::type',
                                   '::integrators::transforms::Korobov<4>::type',
                                   '::integrators::transforms::Sidi<8>::type'
                              ]
        for transform, cpp_transform in zip(transforms_to_test, target_replacements):
            qmc = Qmc(transform=transform)
            self.assertEqual(qmc.options, {})
            self.assertEqual(qmc.integrator_t, """secdecutil::integrators::Qmc<
                                                                integrand_return_t,
                                                                maximal_number_of_integration_variables,
                                                                %s,
                                                                integrand_t
                                                           >""" % cpp_transform)

    #@attr('active')
    def test_multiple_options(self):
        for fitfunction, formatted_fitfunction in ("none", "None"), ("polysingular", "PolySingular"):
            qmc = Qmc(transform="korobov7", fitfunction=fitfunction)
            self.assertEqual(qmc.integrator_t, """secdecutil::integrators::Qmc<
                                                                integrand_return_t,
                                                                maximal_number_of_integration_variables,
                                                                ::integrators::transforms::Korobov<7>::type,
                                                                integrand_t,
                                                                ::integrators::fitfunctions::%s::type
                                                           >""" % formatted_fitfunction)

        qmc = Qmc(transform="korobov7")
        self.assertEqual(qmc.options,{})

        qmc = Qmc(transform="korobov7", epsrel=1e-7, verbosity=1, fitfunction="polysingular",
                  errormode="largest", generatingvectors="cbcpt_cfftw1_6", seed=4592, devices=[1,2,3])
        target_options = dict(
            epsrel='integrator->epsrel = 1e-07;', verbosity='integrator->verbosity = 1;',
            errormode='integrator->errormode = ::integrators::ErrorMode::largest;',
            generatingvectors='integrator->generatingvectors = ::integrators::generatingvectors::cbcpt_cfftw1_6();',
            seed='integrator->randomgenerator.seed(4592);', devices='integrator->devices = { 1, 2, 3 };'
        )
        self.assertEqual(qmc.options,target_options)
