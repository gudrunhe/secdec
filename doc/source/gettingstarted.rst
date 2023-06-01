Getting Started
===============

After installation, you should have a folder `examples` in your main `pySecDec` directory.
Here we describe a few of the examples available in the `examples` directory.
A full list of examples is given in :ref:`list_of_examples`.

.. _a_simple_example:

A Simple Example
----------------

We first show how to compute a simple dimensionally regulated integral:

.. math::

    \int_0^1 \mathrm{d} x \int_0^1 \mathrm{d} y \ (x+y)^{-2+\epsilon}.


To run the example change to the `easy` directory and run the commands::

    $ python3 generate_easy.py
    $ make -C easy
    $ python3 integrate_easy.py

Additional build options are discussed in the :ref:`next section <building_the_cpp_lib>`. This will evaluate and print the result of the integral::

    Numerical Result: + ((1.00000000000000022e+00,0.00000000000000000e+00) +/- (5.65352153979095401e-17,0.00000000000000000e+00))*eps^-1 + ((3.06852819440053548e-01,0.00000000000000000e+00) +/- (1.18502493127591741e-15,0.00000000000000000e+00)) + O(eps)
    Analytic Result: + (1.000000)*eps^-1 + (0.306853) + O(eps)

The file ``generate_easy.py`` defines the integral and calls `pySecDec` to perform the sector decomposition.
When run it produces the directory `easy` which contains the code required to numerically evaluate the integral.
The make command builds this code and produces a library.
The file ``integrate_easy.py`` loads the integral library and evaluates the integral.
The user is encouraged to copy and adapt these files to evaluate their own integrals.

.. note::

    If the user is interested in evaluating a loop integral there are many convenience functions that make this much easier. Please see :ref:`evaluating_a_loop_integral` for more details.


In ``generate_easy.py`` we first import :func:`make_package <pySecDec.make_package>`, a function which can decompose, subtract and expand regulated integrals and write a C++ package to evaluate them.
To define our integral we give it a `name` which will be used as the name of the output directory and C++ namespace.
The `integration_variables` are declared along with a list of the name of the `regulators`.
We must specify a list of the `requested_orders` to which `pySecDec` should expand our integral in each regulator.
Here we specify ``requested_orders = [0]`` which instructs :func:`make_package <pySecDec.make_package>` to expand the integral up to and including :math:`\mathcal{O}(\epsilon)`.
Next, we declare the `polynomials_to_decompose`, here `sympy` syntax should be used.

.. literalinclude:: ../../examples/easy/generate_easy.py
   :language: python

Once the C++ library has been written and built we run ``integrate_easy.py``.
Here the library is loaded using :class:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>`.
Calling the instance of :class:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>` with ``easy_integral()`` numerically evaluates the integral and returns the result.

.. literalinclude:: ../../examples/easy/integrate_easy.py
   :language: python

.. _evaluating_a_loop_integral:

Evaluating a Loop Integral
--------------------------

A simple example of the evaluation of a loop integral with `pySecDec` is `box1L`.
This example computes a one-loop box with one off-shell leg (with off-shellness ``s1``) and one internal massive line (with mass squared ``msq``), it is shown in :numref:`box1L_diagram`.

.. _box1L_diagram:

.. figure:: _static/box1L.*
    :align: center
    :alt: Diagrammatic representation of `box1L`

    Diagrammatic representation of `box1L`

To run the example change to the `box1L` directory and run the commands::

    $ python3 generate_box1L.py
    $ make -C box1L
    $ python3 integrate_box1L.py

This will print the result of the integral evaluated with Mandelstam invariants ``s=4.0``, ``t=-0.75`` and ``s1=1.25``, ``msq=1.0``::

    eps^-2: -0.142868356275422825 - 1.63596224151119965e-6*I +/- ( 0.00118022544307414272 + 0.000210769456586696187*I )
    eps^-1: 0.639405625715768089 + 1.34277036689902802e-6*I +/- ( 0.00650722394065588166 + 0.000971496627153705891*I )
    eps^0 : -0.425514350373418893 + 1.86892487760861536*I +/- ( 0.00706834403694714484 + 0.0186497890361357298*I )

The file ``generate_box1L.py`` defines the loop integral and calls `pySecDec` to perform the sector decomposition. When run it produces the directory `box1L` which contains the code required to numerically evaluate the integral. The make command builds this code and produces a library. The file ``integrate_box1L.py`` loads the integral library and evaluates the integral for a specified numerical point.

The content of the python files is described in detail in the following sections. The user is encouraged to copy and adapt these files to evaluate their own loop integrals.

Defining a Loop Integral
^^^^^^^^^^^^^^^^^^^^^^^^

To explain the input format, let us look at ``generate_box1L.py`` from the one-loop box example. The first line reads

.. code::

    import pySecDec as psd

This line specifies that the module `pySecDec` should be imported with the alias `psd`.
The functions :func:`loop_package <pySecDec.loop_integral.loop_package>` and :mod:`loop_integral <pySecDec.loop_integral>` from this module will be needed shortly.


The following part contains the definition of the loop integral ``li``:

.. code::

    li = psd.LoopIntegralFromGraph(
        # Give adjacency list and indicate whether the propagator
        # connecting the numbered vertices is massive or massless
        # in the first entry of each list item.
        internal_lines = [['m',[1,2]], ['0',[2,3]], ['0',[3,4]], ['0',[4,1]]],
        # List the names of the external momenta and the labels
        # of the vertecies they are attached to.
        external_lines = [['p1',1], ['p2',2], ['p3',3], ['p4',4]],
        # Define the kinematics and the names of the kinematic
        # invariants.
        replacement_rules = [
            ('p1*p1', 's1'),
            ('p2*p2', '0'),
            ('p3*p3', '0'),
            ('p4*p4', '0'),
            ('p3*p2', 't/2'),
            ('p1*p2', 's/2-s1/2'),
            ('p1*p4', 't/2-s1/2'),
            ('p2*p4', 's1/2-t/2-s/2'),
            ('p3*p4', 's/2'),
            ('m**2', 'msq')
        ]
    )

Here the class :class:`LoopIntegralFromGraph <pySecDec.loop_integral.LoopIntegralFromGraph>` is used to Feynman parametrize the loop integral given the adjacency list. Alternatively, the class :class:`LoopIntegralFromPropagators <pySecDec.loop_integral.LoopIntegralFromPropagators>` can be used to construct the Feynman integral given the momentum representation.

The symbols for the kinematic invariants and the masses also need to be given as an ordered list.
The ordering is important as the numerical values assigned to these list elements at the numerical evaluation stage should have the same order.

.. code::

    Mandelstam_symbols = ['s','t','s1']
    mass_symbols = ['msq']


Next, the function :func:`loop_package <pySecDec.loop_integral.loop_package>` is called. It will create a folder called `box1L`.
It performs the algebraic sector decomposition steps and writes a package containing the C++ code for the numerical evaluation.
The argument `requested_orders` specifies the order in the regulator to which the integral should be expanded.
For a complete list of possible options see  :func:`loop_package <pySecDec.loop_integral.loop_package>`.

.. code::

    psd.loop_package(
        name = 'box1L',

        loop_integral = li,

        real_parameters = Mandelstam_symbols + mass_symbols,

        # the highest order of the final epsilon expansion -->
        # change this value to whatever you think is appropriate
        requested_orders = [0],

        # the method to be used for the sector decomposition
        # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
        decomposition_method = 'geometric'
    )

.. _building_the_integration_lib:

Building the integration Library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After running the python script ``generate_box1L.py`` the folder ``box1L`` is created and should contain the following files and subdirectories

.. code::

    Makefile       box1L.hpp           integrate_box1L.cpp  disteval/
    Makefile.conf  box1L.pdf           box1L_data/          pylink/
    README         integral_names.txt  box1L_integral/      src/

In the folder ``box1L``, typing

.. code::

    $ make

will create the static library ``box1L_integral/libbox1L_integral.a`` and the shared library ``box1L_pylink.so`` which can be linked to external programs.
The ``make`` command can also be run in parallel by using the ``-j`` option. The number of threads each instance of ``tform`` uses can be
set via the environment variable `FORMTHREADS`.

.. versionadded:: 1.4
    The environment variable `FORMOPT` sets FORM's code optimization level. If not set, the value that was passed to :func:`make_package <pySecDec.make_package>`
    or :func:`loop_package <pySecDec.loop_integral.loop_package>` is used.

To build the dynamic library ``libbox1L.so`` set ``dynamic`` as build target:

.. code::

    $ make dynamic

.. versionadded:: 1.6
   To build the *disteval* library (which will consist of multiple files in ``box1L/disteval/`` directory) set ``disteval`` as build target:

.. code::

    $ make disteval

Note that *disteval* libraries are designed with a focus on optimization for modern processors by making use of `AVX2 <https://en.wikipedia.org/wiki/AVX2>`_ and `FMA <https://en.wikipedia.org/wiki/FMA_instruction_set>`_ instruction sets.
For CPUs that support these, best performance is achieved by using the newest compiler available on the system (chosen via the ``CXX`` variable), and by enabling the support of AVX2 and FMA (via the ``CXXFLAGS`` variable).
For example:

.. code::

    $ make disteval CXX="g++-12" CXXFLAGS="-mavx2 -mfma"

To build the libraries with NVidia C Compiler (NVCC) for GPU support, type

.. code::

    $ make SECDEC_WITH_CUDA_FLAGS="-arch=sm_XX"

where ``sm_XX`` must be replaced by the target NVidia GPU architechtures; see the `arch option of NVCC <http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/#options-for-steering-gpu-code-generation>`_.
The ``SECDEC_WITH_CUDA_FLAGS`` variable, which enables GPU code compilation, contains flags which are passed to NVCC during code compilation and linking.
Multiple GPU architectures may be specified as described in the `NVCC manual <http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/#options-for-steering-gpu-code-generation>`_, for example
``SECDEC_WITH_CUDA_FLAGS="-gencode arch=compute_XX,code=sm_XX -gencode arch=compute_YY,code=sm_YY"`` where ``XX`` and ``YY`` are the target GPU architectures. The script 
``examples/easy/print-cuda-arch.sh`` can be used to obtain the compute architecture of your current machine.  

To evaluate the integral numerically a program can call one of these libraries.
How to do this interactively or via a python script is explained in the section :ref:`Python Interface <python_interface>`.
Alternatively, a C++ program can be produced as explained in the section :ref:`C++ Interface <cpp_interface>`.

..  _python_interface:

Python Interface
^^^^^^^^^^^^^^^^

To evaluate the integral for a given numerical point we can use ``integrate_box1L.py``.
First it imports the necessary python packages and loads the C++ library.

.. code::

    from pySecDec.integral_interface import IntegralLibrary
    import sympy as sp

    # load c++ library
    box1L = IntegralLibrary('box1L/box1L_pylink.so')

Next, an integrator is configured for the numerical integration. The full list of available integrators and their options is given in :mod:`integral_interface <pySecDec.integral_interface>`.

.. code::

    # choose integrator
    box1L.use_Vegas(flags=2) # ``flags=2``: verbose --> see Cuba manual

If you want to use GPUs, change to the :mod:`CudaQmc<pySecDec.integral_interface.CudaQmc>` integrator. For example, to run on all available GPUs and CPU cores
using the Korobov transform with weight 3, change the above lines to

.. code::

    # choose integrator
    box1L.use_Qmc(transform='Korobov3')


Calling the ``box1L`` library numerically evaluates the integral.
Note that the order of the real parameters must match that specified in ``generate_box1L.py``.
A list of possible settings for the library, in particular details of how to set the contour deformation parameters, is given in :class:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>`.
To change the accuracy settings of the integration, the most important parameters are ``epsrel``, ``epsabs`` and ``maxeval``, which
can be added to the integrator argument list:

.. code::

    # choose integrator
    box1L.use_Vegas(flags=2,epsrel=0.01, epsabs=1e-07, maxeval=1000000)

.. code::

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box1L(real_parameters=[4.0, -0.75, 1.25, 1.0])

In case of a sign check error (sign_check_error), the arguments ``number_of_presamples``, ``deformation_parameters_maximum``, and ``deformation_parameters_minimum`` as described in
:class:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>` can be used to modify the contour.
At this point the string ``str_integral_with_prefactor`` contains the full result of the integral and can be manipulated as required.
In the ``integrate_box1L.py`` an example is shown how to parse the expression with `sympy` and access individual orders of the regulator.

.. note::

   Instead of parsing the result, it can simply be printed with the line ``print(str_integral_with_prefactor)``.

.. code::

    # convert complex numbers from c++ to sympy notation
    str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
    str_prefactor = str_prefactor.replace(',','+I*')
    str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')

    # convert result to sympy expressions
    integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    prefactor = sp.sympify(str_prefactor)
    integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
    integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

An example of how to loop over several kinematic points is shown in the example `integrate_box1L_multiple_points.py`.

Command-line interface with *disteval*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. versionadded:: 1.6

The *disteval* library, once built, can be used directly from the command line via the `pySecDec.disteval` Python module:

.. code::

    $ python3 -m pySecDec.disteval box1L/disteval/box1L.json s=4.0 t=-0.75 s1=1.25 msq=1.0
    [...]
    [
      (
        +eps^-2*(-1.4285714285714279e-01+9.0159338621354360e-18j)
        +eps^-2*(+3.4562234592930473e-17+1.4290950719747641e-17j)*plusminus
        +eps^-1*(+6.3843370937935406e-01+2.5048341561937569e-10j)
        +eps^-1*(+4.4293092326873179e-10+4.6245608965315405e-10j)*plusminus
        +eps^0*(-4.2634981062934296e-01+1.8664974523210687e+00j)
        +eps^0*(+5.5826851229628189e-06+4.8099553795389634e-06j)*plusminus
      )
    ]

Note that the output is a list of expressions; here a list of a single item.
This is because as we shall see in :ref:`evaluating_a_weighted_sum_of_integrals`, a single library can produce multiple resulting expressions.

The general usage of the command-line interface is:

.. code::

    $ python3 -m pySecDec.disteval integrand.json [options] <var>=value ...

The evaluation can be controlled via the provided command-line options:

* ``--epsabs=<number>``: stop if this absolute precision is reached (default: ``1e-10``);
* ``--epsrel=<number>``: stop if this relative precision is reached (default: ``1e-4``);
* ``--timeout=<number>``: stop after at most this many seconds (defaul: ``inf``);
* ``--points=<number>``: use this initial Quasi-Monte-Carlo lattice size (default: ``1e4``);
* ``--presamples=<number>``: use this many points for presampling (default: ``1e4``);
* ``--shifts=<number>``: use this many lattice shifts per integral (default: ``32``);
* ``--coefficients=<path>``: use coefficients from this directory;
* ``--format=<path>``: output the result in this format (``sympy``, ``mathematica``, or ``json``; default: ``sympy``).

This list of options can also be obtained from within the command line by running:

.. code::

    $ python3 -m pySecDec.disteval --help

Python interface with *disteval*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similarly to :ref:`python_interface`, *disteval* libraries can be used from Python, this time via the :class:`DistevalLibrary <pySecDec.integral_interface.DistevalLibrary>` class.
An example of this usage be found in ``integrate_box1L_disteval.py``.
The example starts by importing the necessary packages and loading the library:

.. code::

    from pySecDec.integral_interface import DistevalLibrary
    import sympy as sp

    box1L = DistevalLibrary('box1L/disteval/box1L.json')

Then, calling the ``box1L`` library to perform the evaluation at the given parameter values:

.. code::

    # integrate
    str_result = box1L(parameters={"s": 4.0, "t": -0.75, "s1": 1.25, "msq": 1.0}, verbose=False)

And finally, converting the result to a `sympy` object and printing it:

.. code::

    # convert result to sympy expressions
    result = sp.sympify(str_result)
    value = result[0].subs({"plusminus": 0})
    error = result[0].coeff("plusminus")

    # examples how to access individual orders
    print('Numerical Result')
    print('eps^-2:', value.coeff('eps',-2), '+/- (', error.coeff('eps',-2), ')')
    print('eps^-1:', value.coeff('eps',-1), '+/- (', error.coeff('eps',-1), ')')
    print('eps^0 :', value.coeff('eps',0), '+/- (', error.coeff('eps',0), ')')

..  _cpp_interface:

C++ Interface
^^^^^^^^^^^^^

Usually it is easier to obtain a numerical result using the :ref:`Python Interface <python_interface>`.
However, the library can also be used directly from C++.
Inside the generated `box1L` folder the file ``integrate_box1L.cpp`` demonstrates this.

After the lines parsing the input parameters, an :cpp:class:`secdecutil::Integrator` is constructed and its parameters are set:

.. code-block:: c++

        // Set up Integrator
        secdecutil::integrators::Qmc<
                                    box1L::integrand_return_t,
                                        box1L::maximal_number_of_integration_variables,
                                        integrators::transforms::Korobov<3>::type,
                                        box1L::user_integrand_t
                                    > integrator;
        integrator.verbosity = 1;

The amplitude is constructed via a call to :cpp:func:`name::make_amplitudes` and packed into a :cpp:type:`name::handler_t`.

.. code-block:: c++

        // Construct the amplitudes
        std::vector<box1L::nested_series_t<box1L::sum_t>> unwrapped_amplitudes =
            box1L::make_amplitudes(real_parameters, complex_parameters, "box1L_coefficients", integrator);

        // Pack amplitudes into handler
        box1L::handler_t<box1L::amplitudes_t> amplitudes
        (
            unwrapped_amplitudes,
            integrator.epsrel, integrator.epsabs
            // further optional arguments: maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
        );
        amplitudes.verbose = true;

If desired, the contour deformation can be adjusted via additional arguments to :cpp:type:`name::handler_t`.

.. seealso::

        :numref:`chapter_secdecutil_amplitude` and :numref:`generated_cpp_amplitude_sum_libs` for more detailed information about :cpp:func:`name::make_amplitudes` and :cpp:type:`name::handler_t`.

To numerically integrate the sum of sectors, the :cpp:type:`name::handler_t::evaluate()` function is called:

.. code-block:: c++

        // compute the amplitudes
        const std::vector<box1L::nested_series_t<secdecutil::UncorrelatedDeviation<box1L::integrand_return_t>>> result = amplitudes.evaluate();

The remaining lines print the result:

.. code-block:: c++

        // print the result
        for (unsigned int amp_idx = 0; amp_idx < box1L::number_of_amplitudes; ++amp_idx)
            std::cout << "amplitude" << amp_idx << " = " << result.at(amp_idx) << std::endl;


The C++ program can be built with the command::

    $ make integrate_box1L

A kinematic point must be specified when calling the ``integrate_box1L`` executable, the input format is::

    $ ./integrate_box1L 4.0 -0.75 1.25 1.0 

where the arguments are the ``real_parameters`` values for (``s``, ``t``, ``s1``, ``msq``).
For integrals depending on ``complex_parameters``, their value is specified by a space separated pair of numbers representing the real and imaginary part.

If your integral is higher than seven dimensional, changing the integral transform to :cpp:type:`integrators::transforms::Baker::type` may improve the accuracy of the result. 
For further options of the QMC integrator we refer to :numref:`chapter_cpp_qmc`.

.. _evaluating_a_weighted_sum_of_integrals:

Evaluating a Weighted Sum of Integrals
--------------------------------------

.. versionadded:: 1.5

Let us examine example ``easy_sum``, which demonstrates how two weighted sums of dimensionally regulated integrals can be evaluated.
The example computes the following two weighted sums:

.. math::

    & 2 s\ I_1 + 3 s\ I_2, \\
    & \frac{s}{2 \epsilon}\ I_1 + \frac{s \epsilon}{3}\ I_2,

where

.. math::

    I_1 & = \int_0^1 \mathrm{d} x \int_0^1 \mathrm{d} y \ (x+y)^{-2+\epsilon}, \\
    I_2 & = \int_0^1 \mathrm{d} x \int_0^1 \mathrm{d} y \ (2x+3y)^{-1+\epsilon}.


First, we import the necessary python packages and open the ``if __name__ == "__main__"`` guard, as required by `multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_.

.. code::

    #!/usr/bin/env python3
    import pySecDec as psd

    if __name__ == "__main__":

The common arguments for the integrals are collected in the ``common_args`` dictionary.

.. code::

        common_args = {
            'real_parameters': ['s'],
            'regulators': ['eps'],
            'requested_orders': [0]
        }

Next, the coefficients of the integrals for each weighted sum are specified.
Each coefficient is specified as a string with an arbitrary arithmetic (i.e. rational) expression.
Coefficients can also be specified as instances of the :class:`Coefficient <pySecDec.code_writer.sum_package.Coefficient>` class.
Coefficients can depend on the regulators, the :func:`sum_package <pySecDec.code_writer.sum_package.sum_package>` function will automatically determine the correct orders to which the coefficients and integrals should be expanded in order to obtain the ``requested_orders``.

.. code::

        coefficients = {
            "sum1" : [
                '2*s',       # easy1
                '3*s'        # easy2
            ],
            "sum2" : [
                's/(2*eps)', # easy1
                's*eps/3'    # easy2
            ]
        }


The integrals are specified using the `MakePackage` wrapper function (which has the same arguments as :func:`make_package <pySecDec.code_writer.make_package>`), for loop integrals the `LoopPackage` wrapper may be used (it has the same arguments as :func:`loop_package <pySecDec.loop_integral.loop_package>`).

.. code::

        integrals = [
            psd.MakePackage('easy1',
                integration_variables = ['x','y'],
                polynomials_to_decompose = ['(x+y)^(-2+eps)'],
                **common_args),
            psd.MakePackage('easy2',
                integration_variables = ['x','y'],
                polynomials_to_decompose = ['(2*x+3*y)^(-1+eps)'],
                **common_args)
        ]

Finally, the list of integrals and coefficients are passed to :func:`sum_package <pySecDec.code_writer.sum_package.sum_package>`. This will generate a C++ library which efficiently evaluates both weighted sums of integrals, sharing the results of the integrals between the different sums. 
    
.. code::

        # generate code sum of (int * coeff)
        psd.sum_package('easy_sum', integrals,
            coefficients = coefficients, **common_args)

The generated C++ library can be :ref:`compiled <building_the_cpp_lib>` and called via the :ref:`python <python_interface>` and/or :ref:`C++ <cpp_interface>` interface as described above.

.. _using_expansion_by_regions_generic_integral:

Using Expansion By Regions (Generic Integral)
---------------------------------------------

.. versionadded:: 1.5

The example ``make_regions_ebr`` provides a simple introduction to the expansion by regions functionality within pySecDec.
For a more detailed discussion of expansion by regions see our paper [PSD21]_.

The necessary packages are loaded and the ``if __name__ == "__main__"`` guard is opened.

.. code::

    #!/usr/bin/env python3
    import pySecDec as psd

    if __name__ == "__main__":

Expansion by regions is applied to a generic integral using the :func:`make_regions <pySecDec.make_regions>` function.

.. code::

        regions_generators = psd.make_regions(
            name = 'make_regions_ebr',
            integration_variables = ['x'],
            regulators = ['delta'],
            requested_orders = [0],
            smallness_parameter = 't',
            polynomials_to_decompose = ['(x)**(delta)','(t + x + x**2)**(-1)'],
            expansion_by_regions_order = 0,
            real_parameters = ['t'],
            complex_parameters = [],
            decomposition_method = 'geometric_infinity_no_primary',
            polytope_from_sum_of=[1]
        )

The output of :func:`make_regions <pySecDec.make_regions.make_regions>` can be passed to :func:`sum_package <pySecDec.code_writer.sum_package.sum_package>` in order to generate a C++ library suitable for evaluating the expanded integral.

.. code::

        psd.sum_package(
            'make_regions_ebr',
            regions_generators,
            regulators = ['delta'],
            requested_orders = [0],
            real_parameters = ['t']
        )

The generated C++ library can be :ref:`compiled <building_the_cpp_lib>` and called via the :ref:`python <python_interface>` and/or :ref:`C++ <cpp_interface>` interface as described above.

.. _using_expansion_by_regions_loop_integral:

Using Expansion By Regions (Loop Integral)
------------------------------------------

.. versionadded:: 1.5

The example ``generate_box1L_ebr`` demonstrates how expansion by regions can be applied to loop integrals within pySecDec by applying it to the 1-loop box integral as described in Section 4.2 of [Mis18]_.
For a more detailed discussion of expansion by regions see our paper [PSD21]_.

First, the necessary packages are loaded and the ``if __name__ == "__main__"`` guard is opened.

.. code::

    #!/usr/bin/env python3

    import pySecDec as psd

    # This example is the one-loop box example in Go Mishima's paper arXiv:1812.04373

    if __name__ == "__main__":

The loop integral can be constructed via the convenience functions in :mod:`loop_integral <pySecDec.loop_integral>`, here we use :class:`LoopintegralFromGraph <pySecDec.loop_integral.LoopIntegralFromGraph>`.
Note that ``powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5"]``, here ``n1`` is an extra regulator required to regulate the singularities which appear when expanding this loop integral.
We use the "trick" of introducing only a single regulator divided by different prime numbers for each power, rather than unique regulators for each propagator (though this is also supported by pySecDec). 
Poles in the extra regulator ``n1`` may appear in individual regions but are expected to cancel when all regions are summed.

.. code::

        # here we define the Feynman diagram
        li = psd.LoopIntegralFromGraph(
            internal_lines = [['mt',[3,1]],['mt',[1,2]],['mt',[2,4]],['mt',[4,3]]],
            external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],
            powerlist=["1+n1","1+n1/2","1+n1/3","1+n1/5"],
            regulators=["eps","n1"],
            # renames the parameters to get the same polynomials as in 1812.04373
            Feynman_parameters=["x%i" % i for i in range(1,5)],

        replacement_rules = [
            # note that in those relations all momenta are incoming
            # general relations:
            ('p1*p1', 'm1sq'),
            ('p2*p2', 'm2sq'),
            ('p3*p3', 'm3sq'),
            ('p4*p4', 'm4sq'),
            ('p1*p2', 's/2-(m1sq+m2sq)/2'),
            ('p1*p3', 't/2-(m1sq+m3sq)/2'),
            ('p1*p4', 'u/2-(m1sq+m4sq)/2'),
            ('p2*p3', 'u/2-(m2sq+m3sq)/2'),
            ('p2*p4', 't/2-(m2sq+m4sq)/2'),
            ('p3*p4', 's/2-(m3sq+m4sq)/2'),
            ('u', '(m1sq+m2sq+m3sq+m4sq)-s-t'),
            # relations for our specific case:
            ('mt**2', 'mtsq'),
            ('m1sq',0),
            ('m2sq',0),
            ('m3sq','mHsq'),
            ('m4sq','mHsq'),
            ('mHsq', 0),
        ])

Expansion by regions is applied to a loop integral using the :func:`loop_regions <pySecDec.loop_integral.loop_regions>` function.
We expand around a small mass `mtsq`.

.. code::

        # find the regions
        terms = psd.loop_regions(
            name = "box1L_ebr",
            loop_integral=li,
            smallness_parameter = "mtsq",
            expansion_by_regions_order=0)

The output of :func:`loop_regions <pySecDec.loop_integral.loop_regions>` can be passed to :func:`sum_package <pySecDec.code_writer.sum_package.sum_package>` in order to generate a C++ library suitable for evaluating the expanded integral.

.. code::

        # write the code to sum up the regions
        psd.sum_package(
            "box1L_ebr",
            terms,
            li.regulators,
            requested_orders = [0,0],
            real_parameters = ['s','t','u','mtsq'],
            complex_parameters = [])


The generated C++ library can be :ref:`compiled <building_the_cpp_lib>` and called via the :ref:`python <python_interface>` and/or :ref:`C++ <cpp_interface>` interface as described above.

.. _list_of_examples:

List of Examples
----------------

Here we list the available examples. For more details regarding each example see [PSD17]_, [PSD18]_, [PSD21]_ and [PSD23]_.

+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **easy**:                  | a simple parametric integral, described in :numref:`a_simple_example`                                                          |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **box1L**:                 | a simple 1-loop, 4-point, 4-propagator integral, described in :numref:`evaluating_a_loop_integral`                             |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **triangle2L**:            | a 2-loop, 3-point, 6-propagator diagram, also known as `P126`                                                                  |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **box2L_numerator**:       | a massless planar on-shell 2-loop, 4-point, 7-propagator box with a numerator, either defined as an inverse propagator         |
|                            | ``box2L_invprop.py`` or in terms of contracted Lorentz vectors ``box2L_contracted_tensor.py``                                  |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **pentabox_fin**:          | a 2-loop, 5-point, 8-propagator diagram, evaluated in :math:`6-2 \epsilon` dimensions where it is finite                       |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **triangle3L**:            | a 2-loop, 3-point, 7-propagator integral, demonstrates that the symmetry finder can significantly reduce the number of sectors |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **formfactor4L**:          | a single-scale 4-loop 3-point integral in :math:`6-2 \epsilon` dimensions                                                      |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **bubble6L**:              | a single-scale 6-loop 2-point integral, evaluated at a Euclidean phase-space point                                             |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **elliptic2L_euclidean**:  | an integral known to contain elliptic functions, evaluated at a Euclidean phase-space point                                    |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **elliptic2L_physical**:   | an integral known to contain elliptic functions, evaluated at a physical phase-space point                                     |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **banana_3mass**:          | a 3-loop 2-point integral with three different internal masses known to contain hyperelliptic functions,                       |
|                            | evaluated at a physical phase-space point                                                                                      |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **hyperelliptic**:         | a 2-loop 4-point nonplanar integral known to contain hyperelliptic functions, evaluated at a physical phase-space point        |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **triangle2L_split**:      | a 2-loop, 3-point, 6-propagator integral without a Euclidean region due to special kinematics                                  |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **Nbox2L_split**:          | three 2-loop, 4-point, 5-propagator integrals that need ``split=True`` due to special kinematics                               |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **hypergeo5F4**:           | a general dimensionally regulated parameter integral                                                                           |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **hz2L_nonplanar**:        | a 2-loop, 4-point, 7-propagator integral with internal and external masses                                                     |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **box1L_ebr**:             | uses expansion by regions to expand a 1-loop box with a small internal mass, this integral is also considered in Section 4.2   |
|                            | of [Mis18]_                                                                                                                    |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **bubble1L_ebr**:          | uses expansion by regions to expand a 1-loop, 2-point integral in various limits, demonstrates the use of an additional        |
|                            | regulator as described in [PSD21]_                                                                                             |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **bubble1L_dotted_ebr**:   | uses expansion by regions to expand a 1-loop, 2-point integral, demonstrates the :math:`t` and :math:`z` methods described in  |
|                            | [PSD21]_                                                                                                                       | 
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **bubble2L_largem_ebr**:   | uses expansion by regions to expand a 1-loop, 2-point integral with a large mass                                               |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **bubble2L_smallm_ebr**:   | uses expansion by regions to expand a 1-loop, 2-point integral with a small mass                                               |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **formfactor1L_ebr**:      | uses expansion by regions to compute various 1-loop, 3-point form factor integrals from the literature, demonstrates the use   |
|                            | of ``add_monomial_regulator_power`` to introduce an additional regulator as described in [PSD21]_                              |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **triangle2L_ebr**:        | uses expansion by regions to compute a 2-loop, 3-point integral with a large mass                                              |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **make_regions_ebr**:      | uses expansion by regions to compute a simple generic integral with a small parameter                                          |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **easy_sum**:              | calculates the sum of two integrals with different coefficients, demonstrates the use of ``sum_package``                       |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **yyyy1L**:                | calculates a 1-loop 4-photon helicity amplitude, demonstrates the use of ``sum_package``                                       |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **two_regulators**:        | an integral involving poles in two different regulators.                                                                       |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **userdefined_cpp**:       | a collection of examples demonstrating how to combine polynomials to be decomposed with other user-defined functions           |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **regions**:               | prints a list of the regions obtained by applying expansion by regions to formfactor1L_massless                                |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **region_tools**:          | demonstrates the standalone usage of suggested_extra_regulator_exponent, extra_regulator_constraints and find_regions          |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
