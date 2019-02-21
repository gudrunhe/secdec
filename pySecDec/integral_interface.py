"""
Integral Interface
------------------

An interface to libraries generated by
:func:`pySecDec.code_writer.make_package` or
:func:`pySecDec.loop_integral.loop_package`.

"""

from ctypes import CDLL, c_void_p, c_char_p, c_bool, c_int, c_uint, c_longlong, c_double, c_ulonglong
from threading import Thread
try:
    from Queue import Queue
except ImportError:
    from queue import Queue

# assuming
# enum qmc_transform_t : int
# {
#     no_transform = -1,
#
#     baker = -2,
#
#     korobov1x1 = 1, korobov1x2 = 2, korobov1x3 = 3, korobov1x4 = 4, korobov1x5 = 5, korobov1x6 = 6,
#     korobov2x1 = 7, korobov2x2 = 8, korobov2x3 = 9, korobov2x4 = 10, korobov2x5 = 11, korobov2x6 = 12,
#     korobov3x1 = 13, korobov3x2 = 14, korobov3x3 = 15, korobov3x4 = 16, korobov3x5 = 17, korobov3x6 = 18,
#     korobov4x1 = 19, korobov4x2 = 20, korobov4x3 = 21, korobov4x4 = 22, korobov4x5 = 23, korobov4x6 = 24,
#     korobov5x1 = 25, korobov5x2 = 26, korobov5x3 = 27, korobov5x4 = 28, korobov5x5 = 29, korobov5x6 = 30,
#     korobov6x1 = 31, korobov6x2 = 32, korobov6x3 = 33, korobov6x4 = 34, korobov6x5 = 35, korobov6x6 = 36,
#
#     sidi1 = -11,
#     sidi2 = -12,
#     sidi3 = -13,
#     sidi4 = -14,
#     sidi5 = -15,
#     sidi6 = -16
# };
known_qmc_transforms = dict(
    none = -1,

    baker = -2,

    korobov1x1 = 1, korobov1x2 = 2, korobov1x3 = 3, korobov1x4 = 4, korobov1x5 = 5, korobov1x6 = 6,
    korobov2x1 = 7, korobov2x2 = 8, korobov2x3 = 9, korobov2x4 = 10, korobov2x5 = 11, korobov2x6 = 12,
    korobov3x1 = 13, korobov3x2 = 14, korobov3x3 = 15, korobov3x4 = 16, korobov3x5 = 17, korobov3x6 = 18,
    korobov4x1 = 19, korobov4x2 = 20, korobov4x3 = 21, korobov4x4 = 22, korobov4x5 = 23, korobov4x6 = 24,
    korobov5x1 = 25, korobov5x2 = 26, korobov5x3 = 27, korobov5x4 = 28, korobov5x5 = 29, korobov5x6 = 30,
    korobov6x1 = 31, korobov6x2 = 32, korobov6x3 = 33, korobov6x4 = 34, korobov6x5 = 35, korobov6x6 = 36,

    sidi1 = -11,
    sidi2 = -12,
    sidi3 = -13,
    sidi4 = -14,
    sidi5 = -15,
    sidi6 = -16
)
for i in range(1,7):
    known_qmc_transforms['korobov' + str(i)] = known_qmc_transforms['korobov%ix%i'%(i,i)]

# assuming
# enum qmc_fitfunction_t : int
# {
#     default_fitfunction = 0,
#
#     none = -1,
#     polysingular = 1
# };
known_qmc_fitfunctions = dict(
    default = 0,

    none = -1,
    polysingular = 1
)

# assuming
# enum qmc_generatingvectors_t: int
# {
#     default_generatingvectors = 0,
#
#     cbcpt_dn1_100 = 1,
#     cbcpt_dn2_6 = 2,
#     cbcpt_cfftw1_6 = 3
# };
known_qmc_generatingvectors = dict(
    default = 0,

    cbcpt_dn1_100 = 1,
    cbcpt_dn2_6 = 2,
    cbcpt_cfftw1_6 = 3
)

class CPPIntegrator(object):
    '''
    Abstract base class for integrators to be used with
    an :class:`.IntegralLibrary`.
    This class holds a pointer to the c++ integrator and
    defines the destructor.

    '''
    def __del__(self):
        if hasattr(self, 'c_integrator_ptr'):
            self.c_lib.free_integrator.restype = None
            self.c_lib.free_integrator.argtypes = [c_void_p]
            self.c_lib.free_integrator(self.c_integrator_ptr)

class MultiIntegrator(CPPIntegrator):
    '''
    .. versionadded:: 1.3.1

    Wrapper for the :cpp:class:`secdecutil::MultiIntegrator`.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    :param low_dim_integrator:
        :class:`CPPIntegrator`;
        The integrator to be used if the integrand is lower
        dimensional than `critical_dim`.

    :param high_dim_integrator:
        :class:`CPPIntegrator`;
        The integrator to be used if the integrand has dimension
        `critical_dim` or higher.

    :param critical_dim:
        integer;
        The dimension below which the `low_dimensional_integrator`
        is used.

    Use this class to switch between integrators based on the
    dimension of the integrand when integrating the `integral_ibrary`.
    For example, ":class:`CQuad` for 1D and :class:`Vegas` otherwise"
    is implemented as::

        integral_library.integrator = MultiIntegrator(integral_library,CQuad(integral_library),Vegas(integral_library),2)

    :class:`MultiIntegrator` can be nested to implement multiple
    critical dimensions. To use e.g. :class:`CQuad` for 1D,
    :class:`Cuhre` for 2D and 3D, and :class:`Vegas` otherwise, do::

        integral_library.integrator = MultiIntegrator(integral_library,CQuad(integral_library),MultiIntegrator(integral_library,Cuhre(integral_library),Vegas(integral_library),4),2)

    .. warning::
        The `integral_library` passed to the integrators must be the
        same for all of them. Furthermore, an integrator can only be
        used to integrate the `integral_library` it has beeen
        constructed with.

    .. warning::
        The :class:`MultiIntegrator` cannot be used with :class:`.CudaQmc`.

    '''
    def __init__(self,integral_library,low_dim_integrator,high_dim_integrator,critical_dim):
        self.low_dim_integrator = low_dim_integrator # keep reference to avoid deallocation
        self.high_dim_integrator = high_dim_integrator # keep reference to avoid deallocation
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_MultiIntegrator.restype = c_void_p
        self.c_lib.allocate_MultiIntegrator.argtypes = [c_void_p, c_void_p, c_int]
        self.c_integrator_ptr = self.c_lib.allocate_MultiIntegrator(low_dim_integrator.c_integrator_ptr,high_dim_integrator.c_integrator_ptr,critical_dim)

class CQuad(CPPIntegrator):
    '''
    Wrapper for the cquad integrator defined in the gsl
    library.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    The other options are defined in :numref:`chapter_cpp_cquad`
    and in the gsl manual.

    '''
    def __init__(self,integral_library,epsrel=1e-2,epsabs=1e-7,n=100,verbose=False,zero_border=0.0):
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_gsl_cquad.restype = c_void_p
        self.c_lib.allocate_gsl_cquad.argtypes = [c_double, c_double, c_uint, c_bool, c_double]
        self.c_integrator_ptr = self.c_lib.allocate_gsl_cquad(epsrel,epsabs,n,verbose,zero_border)

class Vegas(CPPIntegrator):
    '''
    Wrapper for the Vegas integrator defined in the cuba
    library.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    The other options are defined in :numref:`chapter_cpp_cuba`
    and in the cuba manual.

    '''
    def __init__(self,integral_library,epsrel=1e-2,epsabs=1e-7,flags=0,seed=0,mineval=0,maxeval=10**6,zero_border=0.0,nstart=1000,nincrease=500,nbatch=1000,real_complex_together=False):
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_cuba_Vegas.restype = c_void_p
        self.c_lib.allocate_cuba_Vegas.argtypes = [c_double, c_double, c_int, c_int, c_longlong, c_longlong, c_double, c_longlong, c_longlong, c_longlong, c_bool]
        self.c_integrator_ptr = self.c_lib.allocate_cuba_Vegas(epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,nstart,nincrease,nbatch,real_complex_together)

class Suave(CPPIntegrator):
    '''
    Wrapper for the Suave integrator defined in the cuba
    library.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    The other options are defined in :numref:`chapter_cpp_cuba`
    and in the cuba manual.

    '''
    def __init__(self,integral_library,epsrel=1e-2,epsabs=1e-7,flags=0,seed=0,mineval=0,maxeval=10**6,zero_border=0.0,nnew=1000,nmin=10,flatness=25.,real_complex_together=False):
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_cuba_Suave.restype = c_void_p
        self.c_lib.allocate_cuba_Suave.argtypes = [c_double, c_double, c_int, c_int, c_longlong, c_longlong, c_double, c_longlong, c_longlong, c_double, c_bool]
        self.c_integrator_ptr = self.c_lib.allocate_cuba_Suave(epsrel,epsabs,flags,seed,mineval,maxeval,zero_border,nnew,nmin,flatness,real_complex_together)

class Divonne(CPPIntegrator):
    '''
    Wrapper for the Divonne integrator defined in the cuba
    library.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    The other options are defined in :numref:`chapter_cpp_cuba`
    and in the cuba manual.

    '''
    def __init__(self, integral_library, epsrel=1e-2, epsabs=1e-7, flags=0, seed=0, mineval=0, maxeval=10**6,zero_border=0.0,
                                         key1=2000, key2=1, key3=1, maxpass=4, border=0., maxchisq=1.,
                                         mindeviation=.15, real_complex_together=False):
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_cuba_Divonne.restype = c_void_p
        self.c_lib.allocate_cuba_Divonne.argtypes = [c_double, c_double, c_int, c_int, c_longlong, c_longlong,
                                                     c_double, c_int, c_int, c_int, c_int, c_double, c_double,
                                                     c_double, c_bool]
        self.c_integrator_ptr = self.c_lib.allocate_cuba_Divonne(epsrel, epsabs, flags, seed, mineval,maxeval,
                                                                 zero_border, key1, key2, key3, maxpass, border,
                                                                 maxchisq, mindeviation, real_complex_together)

class Cuhre(CPPIntegrator):
    '''
    Wrapper for the Cuhre integrator defined in the cuba
    library.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    The other options are defined in :numref:`chapter_cpp_cuba`
    and in the cuba manual.

    '''
    def __init__(self,integral_library,epsrel=1e-2,epsabs=1e-7,flags=0,mineval=0,maxeval=10**6,zero_border=0.0,key=0,real_complex_together=False):
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_cuba_Cuhre.restype = c_void_p
        self.c_lib.allocate_cuba_Cuhre.argtypes = [c_double, c_double, c_int, c_longlong, c_longlong, c_double, c_int, c_bool]
        self.c_integrator_ptr = self.c_lib.allocate_cuba_Cuhre(epsrel,epsabs,flags,mineval,maxeval,zero_border,key,real_complex_together)

class Qmc(CPPIntegrator):
    '''
    Wrapper for the Qmc integrator defined in the integrators
    library.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    :param errormode:
        string;
        The `errormode` parameter of the Qmc, can be
        ``"default"``, ``"all"``, and ``"largest"``.
        ``"default"`` takes the default from the Qmc
        library. See the Qmc docs for details on the
        other settings.

    :param transform:
        string;
        An integral transform related to the parameter
        `P` of the Qmc. The possible choices
        correspond to the integral transforms of the
        underlying Qmc implementation. Possible values
        are, ``"none"``, ``"baker"``, ``sidi#``,
        ``"korobov#"``, and ``korobov#x#`` where
        any ``#`` (the rank of the Korobov/Sidi transform)
        must be an integer between 1 and 6.

    :param fitfunction:
        string;
        An integral transform related to the parameter
        `F` of the Qmc. The possible choices
        correspond to the integral transforms of the
        underlying Qmc implementation. Possible values
        are ``"default"``, ``"none"``, ``"polysingular"``.

    :param generatingvectors:
        string;
        The name of a set of generating vectors.
        The possible choices correspond to the available generating
        vectors of the underlying Qmc implementation. Possible values
        are ``"default"``, ``"cbcpt_dn1_100"``,
        ``"cbcpt_dn2_6"`` and ``"cbcpt_cfftw1_6"``.

    .. seealso::
        The most important options are described in
        :numref:`chapter_cpp_qmc`.

    The other options are defined in the Qmc docs. If
    an argument is set to 0 then the default of the
    underlying Qmc implementation is used.

    '''
    def __init__(self,integral_library,transform,fitfunction='default',generatingvectors='default',epsrel=0.0,epsabs=0.0,maxeval=0,errormode='default',evaluateminn=0,
                      minn=0,minm=0,maxnperpackage=0,maxmperpackage=0,cputhreads=0,cudablocks=0,cudathreadsperblock=0,verbosity=0,seed=0,devices=[]):
        devices_t = c_int * len(devices)
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_integrators_Qmc.restype = c_void_p
        self.c_lib.allocate_integrators_Qmc.argtypes = [
                                                            c_double, # epsrel
                                                            c_double, # epsabs
                                                            c_ulonglong, # maxeval
                                                            c_int, # errormode
                                                            c_ulonglong, # evaluateminn
                                                            c_ulonglong, # minn
                                                            c_ulonglong, # minm
                                                            c_ulonglong, # maxnperpackage
                                                            c_ulonglong, # maxmperpackage
                                                            c_ulonglong, # cputhreads
                                                            c_ulonglong, # cudablocks
                                                            c_ulonglong, # cudathreadsperblock
                                                            c_ulonglong, # verbosity
                                                            c_longlong, # seed
                                                            c_int, # transform_id
                                                            c_int, # fitfunction_id
                                                            c_int # generatingvectors_id
                                                      ]

        # assuming:
        # enum ErrorMode : int
        # {
        #     all = 1,
        #     largest = 2
        # };
        if errormode == 'default':
            errormode_enum = 0
        elif errormode == 'all':
            errormode_enum = 1
        elif errormode == 'largest':
            errormode_enum = 2
        else:
            raise ValueError('Unknown `errormode` "' + str(errormode) + '"')

        self.c_integrator_ptr = self.c_lib.allocate_integrators_Qmc(epsrel,epsabs,maxeval,errormode_enum,evaluateminn,minn,
                                                                    minm,maxnperpackage,maxmperpackage,cputhreads,
                                                                    cudablocks,cudathreadsperblock,verbosity,
                                                                    seed,known_qmc_transforms[str(transform).lower()],
                                                                    known_qmc_fitfunctions[str(fitfunction).lower()],
                                                                    known_qmc_generatingvectors[str(generatingvectors).lower()]
                                                                   )

class CudaQmc(object):
    '''
    Wrapper for the Qmc integrator defined in the integrators
    library for GPU use.

    :param integral_library:
        :class:`IntegralLibrary`;
        The integral to be computed with this integrator.

    :param errormode:
        string;
        The `errormode` parameter of the Qmc, can be
        ``"default"``, ``"all"``, and ``"largest"``.
        ``"default"`` takes the default from the Qmc
        library. See the Qmc docs for details on the
        other settings.

    :param transform:
        string;
        An integral transform related to the parameter
        `P` of the Qmc. The possible choices
        correspond to the integral transforms of the
        underlying Qmc implementation. Possible values
        are, ``"none"``, ``"baker"``, ``sidi#``,
        ``"korobov#"``, and ``korobov#x#`` where
        any ``#`` (the rank of the Korobov/Sidi transform)
        must be an integer between 1 and 6.

    :param fitfunction:
        string;
        An integral transform related to the parameter
        `F` of the Qmc. The possible choices
        correspond to the integral transforms of the
        underlying Qmc implementation. Possible values
        are ``"default"``, ``"none"``, ``"polysingular"``.

    :param generatingvectors:
        string;
        The name of a set of generating vectors.
        The possible choices correspond to the available generating
        vectors of the underlying Qmc implementation. Possible values
        are ``"default"``, ``"cbcpt_dn1_100"``,
        ``"cbcpt_dn2_6"`` and ``"cbcpt_cfftw1_6"``.

    The other options are defined in the Qmc docs. If
    an argument is set to 0 then the default of the
    underlying Qmc implementation is used.

    '''
    def __init__(self,integral_library,transform,fitfunction='default',generatingvectors='default',epsrel=0.0,epsabs=0.0,maxeval=0,errormode='default',evaluateminn=0,
                      minn=0,minm=0,maxnperpackage=0,maxmperpackage=0,cputhreads=0,cudablocks=0,cudathreadsperblock=0,verbosity=0,seed=0,devices=[]):
        devices_t = c_int * len(devices)
        argtypes = [
                        c_double, # epsrel
                        c_double, # epsabs
                        c_ulonglong, # maxeval
                        c_int, # errormode
                        c_ulonglong, # evaluateminn
                        c_ulonglong, # minn
                        c_ulonglong, # minm
                        c_ulonglong, # maxnperpackage
                        c_ulonglong, # maxmperpackage
                        c_ulonglong, # cputhreads
                        c_ulonglong, # cudablocks
                        c_ulonglong, # cudathreadsperblock
                        c_ulonglong, # verbosity
                        c_longlong, # seed
                        c_int, # transform_id
                        c_int, # fitfunction_id
                        c_int, # generatingvectors_id
                        c_ulonglong, # number_of_devices
                        devices_t # devices[]
                   ]
        self.c_lib = integral_library.c_lib
        self.c_lib.allocate_cuda_integrators_Qmc_together.restype = self.c_lib.allocate_cuda_integrators_Qmc_separate.restype = c_void_p
        self.c_lib.allocate_cuda_integrators_Qmc_together.argtypes = self.c_lib.allocate_cuda_integrators_Qmc_separate.argtypes = argtypes

        # assuming:
        # enum ErrorMode : int
        # {
        #     all = 1,
        #     largest = 2
        # };
        if errormode == 'default':
            errormode_enum = 0
        elif errormode == 'all':
            errormode_enum = 1
        elif errormode == 'largest':
            errormode_enum = 2
        else:
            raise ValueError('Unknown `errormode` "' + str(errormode) + '"')

        self.c_integrator_ptr_together = self.c_lib.allocate_cuda_integrators_Qmc_together(
                                                                                               epsrel,epsabs,maxeval,errormode_enum,evaluateminn,minn,
                                                                                               minm,maxnperpackage,maxmperpackage,cputhreads,
                                                                                               cudablocks,cudathreadsperblock,verbosity,
                                                                                               seed,known_qmc_transforms[str(transform).lower()],
                                                                                               known_qmc_fitfunctions[str(fitfunction).lower()],
                                                                                               known_qmc_generatingvectors[str(generatingvectors).lower()],
                                                                                               len(devices),devices_t(*devices)
                                                                                          )
        self.c_integrator_ptr_separate = self.c_lib.allocate_cuda_integrators_Qmc_separate(
                                                                                               epsrel,epsabs,maxeval,errormode_enum,evaluateminn,minn,
                                                                                               minm,maxnperpackage,maxmperpackage,cputhreads,
                                                                                               cudablocks,cudathreadsperblock,verbosity,
                                                                                               seed,known_qmc_transforms[str(transform).lower()],
                                                                                               known_qmc_fitfunctions[str(fitfunction).lower()],
                                                                                               known_qmc_generatingvectors[str(generatingvectors).lower()],
                                                                                               len(devices),devices_t(*devices)
                                                                                          )

    def __del__(self):
        if hasattr(self, 'c_integrator_ptr_together'):
            self.c_lib.free_cuda_together_integrator.restype = None
            self.c_lib.free_cuda_together_integrator.argtypes = [c_void_p]
            self.c_lib.free_cuda_together_integrator(self.c_integrator_ptr_together)

        if hasattr(self, 'c_integrator_ptr_separate'):
            self.c_lib.free_cuda_separate_integrator.restype = None
            self.c_lib.free_cuda_separate_integrator.argtypes = [c_void_p]
            self.c_lib.free_cuda_separate_integrator(self.c_integrator_ptr_separate)

class IntegralLibrary(object):
    r'''
    Interface to a c++ library produced by
    :func:`.make_package` or :func:`.loop_package`.

    :param shared_object_path:
        str;
        The path to the file "<name>_pylink.so"
        that can be built by the command

        .. code::

            $ make pylink

        in the root directory of the c++ library.

    Instances of this class can be called with the
    following arguments:

    :param real_parameters:
        iterable of float;
        The real_parameters of the library.

    :param complex_parameters:
        iterable of complex;
        The complex parameters of the library.

    :param together:
        bool, optional;
        Whether to integrate the sum of all sectors
        or to integrate the sectors separately.
        Default: ``True``.

    :param number_of_presamples:
        unsigned int, optional;
        The number of samples used for the
        contour optimization.
        A larger value here may resolve a sign
        check error (sign_check_error).
        This option is ignored if the integral
        library was created without deformation.
        Default: ``100000``.

    :param deformation_parameters_maximum:
        float, optional;
        The maximal value the deformation parameters
        :math:`\lambda_i` can obtain.
        Lower this value if you get a sign check
        error (sign_check_error).
        If ``number_of_presamples=0``, all
        :math:`\lambda_i` are set to this value.
        This option is ignored if the integral
        library was created without deformation.
        Default: ``1.0``.

    :param deformation_parameters_minimum:
        float, optional;
        The minimal value the deformation parameters
        :math:`\lambda_i` can obtain.
        Lower this value if you get a sign check
        error (sign_check_error).
        If ``number_of_presamples=0``, all
        :math:`\lambda_i` are set to this value.
        This option is ignored if the integral
        library was created without deformation.
        Default: ``1e-5``.

    :param deformation_parameters_decrease_factor:
        float, optional;
        If the sign check with the optimized
        :math:`\lambda_i` fails during the presampling
        stage, all :math:`\lambda_i` are multiplied
        by this value until the sign check passes.
        We recommend to rather change
        ``number_of_presamples``,
        ``deformation_parameters_maximum``,
        and ``deformation_parameters_minimum``
        in case of a sign check error.
        This option is ignored if the integral
        library was created without deformation.
        Default: ``0.9``.

    The call operator returns three strings:
    * The integral without its prefactor
    * The prefactor
    * The integral multiplied by the prefactor

    The integrator can be configured by calling the
    member methods :meth:`.use_Vegas`, :meth:`.use_Suave`,
    :meth:`.use_Divonne`, :meth:`.use_Cuhre`,
    :meth:`.use_CQuad`, and :meth:`.use_Qmc`.
    The available options are listed in the documentation of
    :class:`.Vegas`, :class:`.Suave`, :class:`.Divonne`,
    :class:`.Cuhre`, :class:`.CQuad`, :class:`.Qmc`
    (:class:`.CudaQmc` for GPU version), respectively.
    :class:`CQuad` can only be used for one dimensional
    integrals. A call to :meth:`use_CQuad` configures the
    integrator to use :class:`CQuad` if possible (1D) and the
    previously defined integrator otherwise.
    By default, :class:`CQuad` (1D only) and :class:`Vegas`
    are used with their default arguments.
    For details about the options, refer to the cuba and the
    gsl manual.

    Further information about the library is stored in
    the member variable `info` of type :class:`dict`.

    '''
    def __init__(self, shared_object_path):
        self._cuda = False

        # import c++ library
        c_lib = self.c_lib = CDLL(shared_object_path)


        # set c prototypes
        c_lib.allocate_string.restype = c_void_p
        c_lib.allocate_string.argtypes = None

        c_lib.free_string.restype = None
        c_lib.free_string.argtypes = [c_void_p]

        c_lib.get_integral_info.restype = c_void_p
        c_lib.get_integral_info.argtypes = [c_void_p]

        c_lib.string2charptr.restype = c_char_p
        c_lib.string2charptr.argtypes = [c_void_p]


        # get integral info
        cpp_str_integral_info = c_lib.allocate_string()
        c_lib.get_integral_info(cpp_str_integral_info)
        str_integral_info = c_lib.string2charptr(cpp_str_integral_info)
        c_lib.free_string(cpp_str_integral_info)
        if not isinstance(str_integral_info, str):
            str_integral_info = str_integral_info.decode('ASCII')


        # store the integral info in a dictionary
        integral_info = self.info = dict()
        for line in str_integral_info.split('\n'):
            key, value = line.split('=')
            integral_info[key.strip()] = value.strip(' ,')


        # continue set c prototypes
        self.real_parameter_t = c_double * int(integral_info['number_of_real_parameters'])
        self.complex_parameter_t = c_double * (2*int(integral_info['number_of_complex_parameters'])) # flattened as: ``real(x0), imag(x0), real(x1), imag(x1), ...``

        c_lib.compute_integral.restype = None
        c_lib.compute_integral.argtypes = [
                                               c_void_p, c_void_p, c_void_p, # output strings
                                               c_void_p, # integrator
                                               self.real_parameter_t, # double array
                                               self.complex_parameter_t, # double array as real(x0), imag(x0), real(x1), imag(x1), ...
                                               c_bool, # together
                                               c_uint, # number_of_presamples
                                               c_double, # deformation_parameters_maximum
                                               c_double, # deformation_parameters_minimum
                                               c_double # deformation_parameters_decrease_factor
                                          ]

        # set cuda integrate types if applicable
        try:
            c_lib.cuda_compute_integral.restype = None
            c_lib.cuda_compute_integral.argtypes = [
                                                        c_void_p, c_void_p, c_void_p, # output strings
                                                        c_void_p, # together integrator
                                                        c_void_p, # separate integrator
                                                        self.real_parameter_t, # double array
                                                        self.complex_parameter_t, # double array as real(x0), imag(x0), real(x1), imag(x1), ...
                                                        c_bool, # together
                                                        c_uint, # number_of_presamples
                                                        c_double, # deformation_parameters_maximum
                                                        c_double, # deformation_parameters_minimum
                                                        c_double # deformation_parameters_decrease_factor
                                                   ]
        except AttributeError:
            # c_lib has been compiled without cuda
            pass

        # set the default integrator (CQuad for 1D, Vegas otherwise)
        self.use_Vegas()
        self.use_CQuad()

    def __call__(
                     self, real_parameters=[], complex_parameters=[], together=True,
                     number_of_presamples=100000, deformation_parameters_maximum=1.,
                     deformation_parameters_minimum=1.e-5,
                     deformation_parameters_decrease_factor=0.9
                ):
        # Initialize and launch the underlying c routines in a subprocess
        # to enable KeyboardInterrupt and avoid crashing the primary python
        # interpreter on error.
        queue = Queue()
        integration_thread = Thread(
                                         target=self._call_implementation,
                                         args=(
                                                  queue, real_parameters,
                                                  complex_parameters, together,
                                                  number_of_presamples,
                                                  deformation_parameters_maximum,
                                                  deformation_parameters_minimum,
                                                  deformation_parameters_decrease_factor
                                              )
                                     )
        integration_thread.daemon = True # daemonize worker to have it killed when the main thread is killed
        integration_thread.start()
        while integration_thread.is_alive(): # keep joining worker until it is finished
            integration_thread.join(5) # call `join` with `timeout` to keep the main thread interruptable
        return queue.get()

    def _call_implementation(
                                self, queue, real_parameters, complex_parameters, together,
                                number_of_presamples, deformation_parameters_maximum,
                                deformation_parameters_minimum,
                                deformation_parameters_decrease_factor
                            ):
        # Passed in correct number of parameters?
        assert len(real_parameters) == int(self.info['number_of_real_parameters']), \
            'Passed %i `real_parameters` but %s needs %i.' % (len(real_parameters),self.info['name'],int(self.info['number_of_real_parameters']))
        assert len(complex_parameters) == int(self.info['number_of_complex_parameters']), \
            'Passed %i `complex_parameters` but `%s` needs %i.' % (len(complex_parameters),self.info['name'],int(self.info['number_of_complex_parameters']))

        # set parameter values
        #   - real parameters
        c_real_parameters = self.real_parameter_t(*real_parameters)

        #   - complex parameters
        flattened_complex_parameters = []
        for c in complex_parameters:
            flattened_complex_parameters.append(c.real)
            flattened_complex_parameters.append(c.imag)
        c_complex_parameters = self.complex_parameter_t(*flattened_complex_parameters)

        # allocate c++ strings
        cpp_str_integral_without_prefactor = self.c_lib.allocate_string()
        cpp_str_prefactor = self.c_lib.allocate_string()
        cpp_str_integral_with_prefactor = self.c_lib.allocate_string()

        # call the underlying c routine
        if self._cuda:
            self.c_lib.cuda_compute_integral(
                                                 cpp_str_integral_without_prefactor,
                                                 cpp_str_prefactor, cpp_str_integral_with_prefactor,
                                                 self.integrator.c_integrator_ptr_together,
                                                 self.integrator.c_integrator_ptr_separate, c_real_parameters,
                                                 c_complex_parameters, together,
                                                 number_of_presamples, deformation_parameters_maximum,
                                                 deformation_parameters_minimum,
                                                 deformation_parameters_decrease_factor
                                            )
        else:
            self.c_lib.compute_integral(
                                            cpp_str_integral_without_prefactor,
                                            cpp_str_prefactor, cpp_str_integral_with_prefactor,
                                            self.integrator.c_integrator_ptr, c_real_parameters,
                                            c_complex_parameters, together,
                                            number_of_presamples, deformation_parameters_maximum,
                                            deformation_parameters_minimum,
                                            deformation_parameters_decrease_factor
                                       )

        # convert c++ stings to python strings or bytes (depending on whether we use python2 or python3)
        str_integral_without_prefactor = self.c_lib.string2charptr(cpp_str_integral_without_prefactor)
        str_prefactor = self.c_lib.string2charptr(cpp_str_prefactor)
        str_integral_with_prefactor = self.c_lib.string2charptr(cpp_str_integral_with_prefactor)

        # free allocated c++ strings
        self.c_lib.free_string(cpp_str_integral_without_prefactor)
        self.c_lib.free_string(cpp_str_prefactor)
        self.c_lib.free_string(cpp_str_integral_with_prefactor)

        # python 2/3 compatibility: make sure the strings read from c++ have type "str" with ASCII encoding
        if not isinstance(str_integral_without_prefactor, str):
            str_integral_without_prefactor = str_integral_without_prefactor.decode('ASCII')
        if not isinstance(str_prefactor, str):
            str_prefactor = str_prefactor.decode('ASCII')
        if not isinstance(str_integral_with_prefactor, str):
            str_integral_with_prefactor = str_integral_with_prefactor.decode('ASCII')

        queue.put( (str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor) )

    def use_Vegas(self, *args, **kwargs):
        self._cuda = False
        self.high_dimensional_integrator = self.integrator = Vegas(self,*args,**kwargs)

    def use_Suave(self, *args, **kwargs):
        self._cuda = False
        self.high_dimensional_integrator = self.integrator = Suave(self,*args,**kwargs)

    def use_Divonne(self, *args, **kwargs):
        self._cuda = False
        self.high_dimensional_integrator = self.integrator = Divonne(self,*args,**kwargs)

    def use_Cuhre(self, *args, **kwargs):
        self.high_dimensional_integrator = self.integrator = Cuhre(self,*args,**kwargs)

    def use_CQuad(self, *args, **kwargs):
        if self._cuda:
            raise RuntimeError('Cannot use `CQuad` together with `CudaQmc`.')
        self.cquad = CQuad(self, *args, **kwargs)
        self.integrator = MultiIntegrator(self,self.cquad,self.high_dimensional_integrator,2)

    def use_Qmc(self, *args, **kwargs):
        if hasattr(self.c_lib,'allocate_integrators_Qmc'):
            self._cuda = False
            self.integrator = Qmc(self,*args,**kwargs)
        else:
            self._cuda = True
            self.integrator = CudaQmc(self,*args,**kwargs)
