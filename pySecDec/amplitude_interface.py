"""
Amplitude Interface
-------------------

An interface to libraries generated by
:func:`pySecDec.code_writer.sum_package`

"""

# TODO: combine with module "integral_interface" into "interface.integral" and "interface.amplitude"

import re
from ctypes import CDLL, c_void_p, c_char_p, c_bool, c_int, c_uint, c_longlong, c_double, c_ulonglong
from threading import Thread
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
import os.path
import collections
import pySecDec

class Configuration(object):
    """
    Class for producing the integrator configuration file "config_<name>.hpp".

    :param main_integral_name:
        string;
        The path to the main integral folder generated by
        :func:`pySecDec.code_writer.sum_package`.

    """
#    It reads in the defined integrals from a file and then only accepts
#    integrators defined for these integrals, otherwise gives an error.
#    The user can create the integrators (i.e. :class:`Vegas`) and set them
#    for an integral using the :func:`Configuration.set_integrator` function,
#    then a call to the :func:`Configuration.configure` function writes the
#    configuration file.
#
#    Example:
#
#    Let's say the integral structure is:
#
#    integral1
#
#    * subintegral1
#    * subintegral2
#
#    To set integrator ``int1=CubaIntegrator("Vegas", epsrel=0.001)`` for every integral under ``integral1``,
#    but ``int2=CubaIntegrator("Cuhre")`` for ``subintegral2``, and then write the
#    config_name.hpp to make pySecDec use these integrators, use:
#
#    .. code-block:: python
#
#        int1 = CubaIntegrator("Vegas", epsrel=0.001)
#        int2 = CubaIntegrator("Cuhre")
#        contourdef = ContourDeformation()
#        conf = Configuration("/path/to/integral/folder")
#        conf.set_options(int1, contourdef)
#        conf.set_options(int2, contourdef, "subintegral2")
#        conf.configure()

    # this is a template for the config_name.hpp file
    config_hpp_template = '''#ifndef %(name)s_config_%(name)s_hpp_included
#define %(name)s_config_%(name)s_hpp_included

/*
 * Set contour-deformation and integrator properties. These can be set per
 * integral by defining an Options struct in the appropriate namespace
 * (see below). Integranls for which no special Options struct is defined
 * use the Options struct in the parent namespace of all the integrands.
 */

%(integrals_config)s

#endif
'''

    # this is a template for the options struct for one integral
    integral_configuration_template = '''// This struct contains the options for the integrals of the enclosing namespace.
struct Options
{
    // integrator options
    // whether the integrator to be used supports CUDA (out of the builtin integrators, only the Qmc does)
    constexpr static bool cuda_compliant_integrator = %(cuda_compliant_integrator)s;

    // integral settings
    template<typename integrand_t, unsigned long long int maximal_number_of_integration_variables>
    struct Integral
    {
        // the integrator (Cuhre, Divonne, Qmc, Suave, Vegas)
        using integrator_t = %(integrator_t)s;

        // instantiation and options of the integrator
        static std::shared_ptr<integrator_t> configure_integrator()
        {
            std::shared_ptr<integrator_t> integrator = std::make_shared<integrator_t>();
            %(integrator_options)s
            return integrator;
        }

        // the kind of integral (CubaIntegral, QmcIntegral, <CustomIntegral>)
        using integral_t = %(integral_t)s;
    };

    // contour deformation parameters
    struct ContourDeformation
    {
        constexpr static unsigned number_of_presamples = %(number_of_presamples)s;
        constexpr static real_t deformation_parameters_maximum = %(deformation_parameters_maximum)s;
        constexpr static real_t deformation_parameters_minimum = %(deformation_parameters_minimum)s;
        constexpr static real_t deformation_parameters_decrease_factor = %(deformation_parameters_decrease_factor)s;
    };
};'''

    def __init__(self, main_integral_name):
        self.integrator_names_file = os.path.join(main_integral_name,"integral_names.txt")

        with open(self.integrator_names_file) as filehandle:
            # read in the integral names
            integral_names = filehandle.read().split()
        assert len(integral_names) > 1, "integral_names.txt must have the name of the main integral and of at least one subintegral"
        self.main_integral_name = integral_names[0]
        self.sub_integral_names = integral_names[1:]

        # for each integral, this dict contains a dict with the integrator and contour deformation info
        # i.e self.integral_options = {"main_integral": {"integrator": integr, "contourdefpar": contdef},
        #       "sub_integral": {"integrator": integr2, "contourdefpar": contdef2}}
        self.integral_options = {name: {} for name in integral_names}

        # this is the file where the configuration will be written
        self.config_filename = os.path.join(main_integral_name, "config_"+self.main_integral_name+".hpp")

    def set_options(self, integrator, contour_deformation, sub_integral_name=None):
        """
        Set the integrator and contour deformation parameters for the integral.

        The options for the main integral apply to all subintegrals unless
        specific options are set for a subintegral. In case options are
        defined for both, the main integral and a subintegral,
        then the subintegral-specific options are used.

        :param integrator:
            :class:`.CubaIntegrator` or :class:`.Qmc`;
            The integrator to be set for the integral.

        :param contour_deformation:
            :class:`pySecDec.amplitude_interface.ContourDeformation`;
            Contains the contour deformation parameters

        :param sub_integral_name:
            string or None;
            Name of the integral to apply the integrator. If ``None``, then
            the settings are applied at the main integral level; i.e. for
            all integrals which do not have specific options defined.

        """

        # when no integral name is given, set it for all the integrals (i.e. the main integral)
        if sub_integral_name is None:
            integral_name = self.main_integral_name
        else:
            integral_name = sub_integral_name

        if integral_name in self.integral_options:
            # the integral name exists, so add the options to the ingeral_options dict
            self.integral_options[integral_name]["integrator"] = integrator
            self.integral_options[integral_name]["contourdefpar"] = contour_deformation
        else:
            # the integral name is not known
            raise NameError("Integral \""+integral_name+"\" does not exist")

    def configure(self):
        """
        Write the config_<name>.hpp file. A c++-options struct is generated
        for each integral that was registered in a call to :meth:`.set_options`.

        """

        # start producing the output text
        options_text = "namespace %s\n{" % self.main_integral_name

        # add the main integral settings if they are set
        main_integral_options = self.integral_options[self.main_integral_name]
        if main_integral_options:
            # the integral has options defined, so add them to the file contents string
            options_text += "\n"+" "*4 + self._make_cpp_options_struct(main_integral_options).replace("\n","\n"+" "*4)+"\n"

        # add the subintegral settings if they are set
        for integral_name in self.sub_integral_names:
            subintegral_options = self.integral_options[integral_name]
            if subintegral_options:
                # the integral has options defined, so add them to the file contents string
                options_text += "\n" + " "*4 + "namespace %s\n    {" % integral_name
                options_text += "\n" + " "*8 + self._make_cpp_options_struct(subintegral_options).replace("\n","\n"+" "*8)
                options_text += "\n" + " "*4 + "}"+"\n"

        options_text += "}"

        # removd whitespace on empty lines
        options_text = re.sub(r"($|\n) +(^|\n)", r"\g<1>\g<2>", options_text)

        # write the generated settings string in the file
        with open(self.config_filename, mode="w") as file:
            file.write(Configuration.config_hpp_template % {"name": self.main_integral_name, "integrals_config": options_text})

    def _make_cpp_options_struct(self, integral_info):
#        """
#        This returns an options string for one integral.
#
#        :param integral_info:
#            dictionary;
#            Contains information about the integral. It should have entries
#            `"integrator"` and `"contourdefpar"` to set the corresponding
#            options.
#
#        Example:
#
#        The following code:
#
#        .. code-block:: python
#
#            conf = Configuration("integralname")
#            string = conf.get_options_struct({"integrator": Qmc(transform="korobov5x5"), "contourdefpar": ContourDeformation()})
#            print(string)
#
#        produces the following output::
#
#            // This struct contains the options for the integrals of the enclosing namespace.
#            struct Options
#            {
#               // integrator options
#               // whether the integrator to be used supports CUDA (out of the builtin integrators, only the Qmc does)
#               constexpr static bool cuda_compliant_integrator = true;
#
#               // integral settings
#               template<typename integrand_t, unsigned long long int maximal_number_of_integration_variables>
#               struct Integral
#               {
#                   // the integrator (Cuhre, Divonne, Qmc, Suave, Vegas)
#                   using integrator_t = secdecutil::integrators::Qmc<
#                                                                           integrand_return_t,
#                                                                           maximal_number_of_integration_variables,
#                                                                           ::integrators::transforms::Korobov<5, 5>::type,
#                                                                           integrand_t
#                                                                   >;
#
#                   // instantiation and options of the integrator
#                   static std::shared_ptr<integrator_t> configure_integrator()
#                   {
#                       std::shared_ptr<integrator_t> integrator = std::make_shared<integrator_t>();
#                       integrator->transform = korobov5x5;
#                       return integrator;
#                   }
#
#                   // the kind of integral (CubaIntegral, QmcIntegral, <CustomIntegral>)
#                   using integral_t = secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,integrator_t,integrand_t>;
#               };
#
#               // contour deformation parameters
#               struct ContourDeformation
#               {
#                   constexpr static unsigned number_of_presamples = 100000;
#                   constexpr static real_t deformation_parameters_maximum = 1.0;
#                   constexpr static real_t deformation_parameters_minimum = 1e-05;
#                   constexpr static real_t deformation_parameters_decrease_factor = 0.9;
#               };
#           };
#
#       """

        integrator = integral_info.get("integrator", False)
        contourdefpar = integral_info.get("contourdefpar", False)
        assert integrator and contourdefpar, "To define an Options struct, both the integrator"\
                                            " and the contour deformation parameters have to be set"

        integrator_options = "// no integrator options set"
        if integrator.options:
            options_to_insert = [optiontext for option,optiontext in integrator.options.items()]
            integrator_options = ("\n"+" "*4*3).join(options_to_insert) # add the proper number of spaces

        replacements = {"cuda_compliant_integrator": str(integrator.cuda_compliant_integrator).lower(),
                        "integrator_t": integrator.integrator_t, "integral_t": integrator.integral_t,
                        "integrator_options": integrator_options,
                        "number_of_presamples": contourdefpar.number_of_presamples,
                        "deformation_parameters_maximum": contourdefpar.deformation_parameters_maximum,
                        "deformation_parameters_minimum": contourdefpar.deformation_parameters_minimum,
                        "deformation_parameters_decrease_factor": contourdefpar.deformation_parameters_decrease_factor
                        }

        return Configuration.integral_configuration_template % replacements

class ContourDeformation(object):
    """
    A class containing the contour deformation parameters. For the meaning
    of the contour deformation parameters, see
    :class:`pySecDec.integral_interface.IntegralLibrary`.

    """
    def __init__(self, number_of_presamples=100000, deformation_parameters_maximum=1.0,
                                deformation_parameters_minimum=1e-5, deformation_parameters_decrease_factor=0.9):
        self.number_of_presamples = number_of_presamples
        self.deformation_parameters_maximum = deformation_parameters_maximum
        self.deformation_parameters_minimum = deformation_parameters_minimum
        self.deformation_parameters_decrease_factor = deformation_parameters_decrease_factor

class CubaIntegrator(object):
    """
    The Cuba integrators: Vegas, Suave, Divonne, and Cuhre.
    For details about the integrator options, see
    :class:`pySecDec.integral_interface.Vegas`,
    :class:`pySecDec.integral_interface.Suave`,
    :class:`pySecDec.integral_interface.Divonne`, and
    :class:`pySecDec.integral_interface.Cuhre`.

    :param integrator_name:
        string;
        The name of the specific integrator. Can be one of the following:
        "Vegas", "Suave", "Divonne", "Cuhre".

    """
    cuda_compliant_integrator = False
    integral_t = "secdecutil::amplitude::CubaIntegral<integrand_return_t,real_t,integrator_t,integrand_t>"

    def __init__(self, integrator_name, **kwargs):
        possible_integrators = "Vegas", "Suave", "Divonne", "Cuhre"
        assert integrator_name in possible_integrators, \
            "The integrator_name must be one of the following: "+", ".join(possible_integrators)+"."
        self.options = {key: "integrator->%s = %s;" % (key,value) for key, value in kwargs.items()}
        self.Integrator_name = integrator_name
        self.integrator_t = "secdecutil::cuba::" + integrator_name + "<integrand_return_t>"

class Qmc(object):
    """
    Qmc integrator. See :class:`pySecDec.integral_interface.Qmc` for
    details about the options.

    """
    cuda_compliant_integrator = True
    integral_t = "secdecutil::amplitude::QmcIntegral<integrand_return_t,real_t,integrator_t,integrand_t>"

    def __init__(self, **kwargs):
        # construct integrator_t
        try:
            transform = kwargs.pop("transform")
        except LookupError:
            raise TypeError("Missing keyword argument `transform`")
        fitfunction = kwargs.pop("fitfunction", None)
        self.integrator_t = self._construct_integrator_t(transform, fitfunction)

        # process other options
        options = self.options = {key: "integrator->%s = %s;" % (key,value) for key, value in kwargs.items()}

        # options which need special treatment
        if "errormode" in kwargs:
            options["errormode"] = "integrator->%s = %s;" % ("errormode", "::integrators::ErrorMode::" + kwargs["errormode"])
        if "generatingvectors" in kwargs:
            options["generatingvectors"] = "integrator->%s = %s;" % \
                ("generatingvectors","::integrators::generatingvectors::" + kwargs["generatingvectors"] + "()")
        if "seed" in kwargs:
            options["seed"] = "integrator->randomgenerator.seed(%i);" % kwargs["seed"]
        if "devices" in kwargs:
            options["devices"] = "integrator->devices = { %s };" % ", ".join(map(str,kwargs["devices"]))

    def _construct_integrator_t(self, transform, fitfunction):
        integrator_t_base = """secdecutil::integrators::Qmc<
                                                                integrand_return_t,
                                                                maximal_number_of_integration_variables,
                                                                %(transform)s,
                                                                integrand_t%(fitfunction)s
                                                           >"""
        # specify fitfunction
        if fitfunction is None:
            cpp_fitfunction = ""
        else:
            for known_fitfunction in "None", "PolySingular": # convert to the proper case, so that user input of "none" and "polysingular" works
                fitfunction = re.sub(known_fitfunction, known_fitfunction, fitfunction, flags=re.IGNORECASE)
            cpp_fitfunction = ",\n" + 16*4*" " + "::integrators::fitfunctions::" + fitfunction + "::type"

        # specify transform
        known_qmc_transforms = {
                                    "none":"::integrators::transforms::None::type",
                                    "baker":"::integrators::transforms::Baker::type",
                                    r"korobov(\d+)x(\d+)":r"::integrators::transforms::Korobov<\g<1>,\g<2>>::type",
                                    r"korobov(\d+)":r"::integrators::transforms::Korobov<\g<1>>::type",
                                    r"sidi(\d+)":r"::integrators::transforms::Sidi<\g<1>>::type"
                               }
        for transform_try in known_qmc_transforms:
            match = re.match(transform_try+r"\Z", transform, flags=re.IGNORECASE)
            if match:
                cpp_transform = match.expand(known_qmc_transforms[transform_try])
                break
        else:
            raise ValueError('Unknown `transform` "' + str(transform) + '"')

        return integrator_t_base % {"transform":cpp_transform, "fitfunction":cpp_fitfunction}
