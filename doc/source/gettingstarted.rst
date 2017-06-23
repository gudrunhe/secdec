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

    $ python generate_easy.py
    $ make -C easy
    $ python integrate_easy.py

This will evaluate and print the result of the integral::

    Numerical Result: + (1.00015897181235158e+00 +/- 4.03392522752491021e-03)*eps^-1 + (3.06903035514056399e-01 +/- 2.82319349818329918e-03) + O(eps)
    Analytic Result: + (1.000000)*eps^-1 + (0.306853) + O(eps)

The file ``generate_easy.py`` defines the integral and calls `pySecDec` to perform the sector decomposition.
When run it produces the directory `easy` which contains the code required to numerically evaluate the integral.
The make command builds this code and produces a library.
The file ``integrate_easy.py`` loads the integral library and evaluates the integral.
The user is encouraged to copy and adapt these files to evaluate their own integrals.

.. note::

    If the user is interested in evaluating a loop integral there are many convenience functions that make this much easier. Please see :ref:`evaluating_a_loop_integral` for more details.


In ``generate_easy.py`` we first import :func:`make_package <pySecDec.code_writer.make_package>`, a function which can decompose, subtract and expand regulated integrals and write a C++ package to evaluate them.
To define our integral we give it a `name` which will be used as the name of the output directory and C++ namespace.
The `integration_variables` are declared along with a list of the name of the `regulators`.
We must specify a list of the `requested_orders` to which `pySecDec` should expand our integral in each regulator.
Here we specify ``requested_orders = [0]`` which instructs :func:`make_package <pySecDec.code_writer.make_package>` to expand the integral up to and including :math:`\mathcal{O}(\epsilon)`.
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

    $ python generate_box1L.py
    $ make -C box1L
    $ python integrate_box1L.py

This will print the result of the integral evaluated with Mandelstam invariants ``s=4.0``, ``t=-0.75`` and ``s1=1.25``, ``msq=1.0``::

    leading pole: -0.142868356275422825 - 1.63596224151119965e-6*I +/- ( 0.00118022544307414272 + 0.000210769456586696187*I )
    subleading pole: 0.639405625715768089 + 1.34277036689902802e-6*I +/- ( 0.00650722394065588166 + 0.000971496627153705891*I )
    finite part: -0.425514350373418893 + 1.86892487760861536*I +/- ( 0.00706834403694714484 + 0.0186497890361357298*I )

The file ``generate_box1L.py`` defines the loop integral and calls `pySecDec` to perform the sector decomposition. When run it produces the directory `box1L` which contains the code required to numerically evaluate the integral. The make command builds this code and produces a library. The file ``integrate_box1L.py`` loads the integral library and evaluates the integral for a specified numerical point.

The content of the python files is described in detail in the following sections. The user is encouraged to copy and adapt these files to evaluate their own loop integrals.

Defining a Loop Integral
^^^^^^^^^^^^^^^^^^^^^^^^

To explain the input format, let us look at ``generate_box1L.py`` from the one-loop box example. The first two lines read

.. code::

    import pySecDec as psd
    from pySecDec.loop_integral import loop_package

They say that the module `pySecDec` should be imported with the alias `psd`, and that the
function :func:`loop_package <pySecDec.loop_integral.loop_package>` from the module :mod:`loop_integral <pySecDec.loop_integral>` is needed.


The following part contains the definition of the loop integral ``li``:

.. code::

    li = psd.loop_integral.LoopIntegralFromGraph(
    # give adjacency list and indicate whether the propagator connecting the numbered vertices is massive or massless in the first entry of each list item.
    internal_lines = [['m',[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
    # contains the names of the external momenta and the label of the vertex they are attached to
    external_lines = [['p1',1],['p2',2],['p3',3],['p4',4]],

    # define the kinematics and the names for the kinematic invariants
    replacement_rules = [
                            ('p1*p1', 's1'),
                            ('p2*p2', 0),
                            ('p3*p3', 0),
                            ('p4*p4', 0),
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
The argument `requested_order` specifies the order in the regulator to which the integral should be expanded.
For a complete list of possible options see  :func:`loop_package <pySecDec.loop_integral.loop_package>`.

.. code::

    loop_package(

    name = 'box1L',

    loop_integral = li,

    real_parameters = Mandelstam_symbols + mass_symbols,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,

    # the optimization level to use in FORM (can be 0, 1, 2, 3)
    form_optimization_level = 2,

    # the WorkSpace parameter for FORM
    form_work_space = '100M',

    # the method to be used for the sector decomposition
    # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    decomposition_method = 'iterative',
    # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
    # $PATH, you can set the path to the 'normaliz' command-line
    # executable here
    #normaliz_executable='/path/to/normaliz',

    )

Building the C++ Library
^^^^^^^^^^^^^^^^^^^^^^^^

After running the python script `generate_box1L.py` the folder `box1L` is created and should contain the following files and subdirectories

.. code::

    Makefile    Makefile.conf    README    box1L.hpp    codegen    integrate_box1L.cpp    pylink    src

in the folder `box1L`, typing

.. code::

    $ make

will create the libraries ``libbox1L.a`` and ``box1L_pylink.so`` which can be linked to an external program calling these integrals.
The ``make`` command can also be run in parallel by using the ``-j`` option.

To evaluate the integral numerically a program can call one of these libraries.
How to do this interactively or via a python script is explained in the section :ref:`Python Interface <python_interface>`.
Alternatively, a C++ program can be produced as explained in the section :ref:`C++ Interface <cpp_interface>`.

..  _python_interface:

Python Interface (basic)
^^^^^^^^^^^^^^^^^^^^^^^^

To evaluate the integral for a given numerical point we can use ``integrate_box1L.py``.
First it imports the necessary python packages and loads the C++ library.

.. code::

    from __future__ import print_function
    from pySecDec.integral_interface import IntegralLibrary
    import sympy as sp

    # load c++ library
    box1L = IntegralLibrary('box1L/box1L_pylink.so')

Next, an integrator is configured for the numerical integration. The full list of available integrators and their options is given in :mod:`integral_interface <pySecDec.integral_interface>`.

.. code::

    # choose integrator
    box.use_Vegas(flags=2) # ``flags=2``: verbose --> see Cuba manual

Calling the ``box`` library numerically evaluates the integral.
Note that the order of the real parameters must match that specified in ``generate_box1L.py``.
A list of possible settings for the library, in particular details of how to set the contour deformation parameters, is given in :class:`IntegralLibrary <pySecDec.integral_interface.IntegralLibrary>`.

.. code::

    # integrate
    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = box1L(real_parameters=[4.0, -0.75, 1.25, 1.0])

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
    print('leading pole:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
    print('subleading pole:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
    print('finite part:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

An example of how to loop over several kinematic points is shown in the example `multiple_kinematic_points.py`.

..  _cpp_interface:

C++ Interface (advanced)
^^^^^^^^^^^^^^^^^^^^^^^^

Usually it is easier to obtain a numerical result using the :ref:`Python Interface <python_interface>`.
However, the library can also be used directly from C++.
Inside the generated `box1L` folder the file ``integrate_box1L.cpp`` demonstrates this.

The function ``print_integral_info`` shows how to access the important variables of the integral library.

In the ``main`` function a kinematic point must be specified by setting the ``real_parameters`` variable, for example::

    int main()
    {
    //  User Specified Phase-space point
        const std::vector<box1L::real_t> real_parameters = {4.0, -0.75, 1.25, 1.0}; // EDIT: kinematic point specified here
        const std::vector<box1L::complex_t> complex_parameters = {  };

The :cpp:func:`name::make_integrands` function returns an :cpp:class:`secdecutil::IntegrandContainer` for each sector and regulator order::

    //  Generate the integrands (optimization of the contour if applicable)
        const std::vector<box1L::nested_series_t<box1L::integrand_t>> sector_integrands = box1L::make_integrands(real_parameters, complex_parameters);

The sectors can be added before integration::

    //  Add integrands of sectors (together flag)
        const box1L::nested_series_t<box1L::integrand_t> all_sectors = std::accumulate(++sector_integrands.begin(), sector_integrands.end(), *sector_integrands.begin() );

An :cpp:class:`secdecutil::Integrator` is constructed and its parameters are set::

    //  Integrate
        secdecutil::cuba::Vegas<box1L::integrand_return_t> integrator;
        integrator.flags = 2; // verbose output --> see cuba manual

To numerically integrate the functions the :cpp:func:`secdecutil::Integrator::integrate` function is applied to each :cpp:class:`secdecutil::IntegrandContainer` using :cpp:func:`secdecutil::deep_apply`::

    const box1L::nested_series_t<secdecutil::UncorrelatedDeviation<box1L::integrand_return_t>> result_all = secdecutil::deep_apply( all_sectors, integrator.integrate );


The remaining lines print the result::

        std::cout << "------------" << std::endl << std::endl;

        std::cout << "-- integral info -- " << std::endl;
        print_integral_info();
        std::cout << std::endl;

        std::cout << "-- integral without prefactor -- " << std::endl;
        std::cout << result_all << std::endl << std::endl;

        std::cout << "-- prefactor -- " << std::endl;
        const box1L::nested_series_t<box1L::integrand_return_t> prefactor = box1L::prefactor(real_parameters, complex_parameters);
        std::cout << prefactor << std::endl << std::endl;

        std::cout << "-- full result (prefactor*integral) -- " << std::endl;
        std::cout << prefactor*result_all << std::endl;
        return 0;
    }

After editing the ``real_parameters`` as described above the C++ program can be built and executed with the commands

.. code::

    $ make integrate_box1L
    $ ./integrate_box1L

.. _list_of_examples:

List of Examples
----------------

Here we list the available examples. For more details regarding each example see [PSD17]_.

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
| **triangle3L**:            | a 2-loop, 3-point, 7-propagator integral, demonstrates that the symmetry finder can significantly reduce the number of sectors |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **elliptic2L_euclidean**:  | an integral known to contain elliptic functions, evaluated at a Euclidean phase-space point                                    |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **elliptic2L_physical**:   | an integral known to contain elliptic functions, evaluated at a physical phase-space point                                     |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **triangle2L_split**:      | a 2-loop, 3-point, 6-propagator integral without a Euclidean region due to special kinematics                                  |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **hypergeo5F4**:           | a general dimensionally regulated parameter integral                                                                           |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **4photon1L_amplitude**:   | calculation of the 4-photon amplitude, showing how to use `pySecDec` as an integral library in a larger context                |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **two_regulators**:        | an integral involving poles in two different regulators.                                                                       |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
| **userdefined_cpp**:       | a collection of examples demonstrating how to combine polynomials to be decomposed with other user-defined functions           |
+----------------------------+--------------------------------------------------------------------------------------------------------------------------------+
