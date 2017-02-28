Getting started
===============

After installation, you should have a folder `examples` in your main pySecDec directory.
It contains various examples, the easiest one being a one-loop box,
`box1L.py`. It also contains some two-loop examples: `triangle2L.py`, 
`box2L.py`, `elliptic_I1.py`, and examples for parametric functions
not related to loop integrals: `Hypergeo5F4.py`
calculates Hypergeomatric functions, which can have (regulated) poles at both zero
and one), `two_regulators.py` contains an example involving poles in two
different regulators. More complex examples are the calcuation of the
4-photon amplitude, which shows how to use pySecDec as an integral
library in a larger context, and the `userdefined_cpp` example which
shows how the user can combine functions to be decomposed with other, user-defined functions. 


User input
----------

To explain the input format, let us look at the one-loop box example. The first two lines read

.. code::

  import pySecDec as psd
  from pySecDec.loop_integral import loop_package

They say that the module `pySecDec` should be imported with the alias `psd`, and that the 
function :func:`loop_package <pySecDec.loop_integral.loop_package>` from the class :class:`LoopIntegral <pySecDec.loop_integral.LoopIntegral>` is needed. 


The following part contains the definition of the loop integral li:

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

The symbols for the kinematic invariants and the masses also need to be given as an ordered list. 
The ordering is important as the numerical values assigned to these list elements at the numerical evaluation stage should have the same order.

.. code::

 Mandelstam_symbols = ['s','t','s1']
 mass_symbols = ['msq']
 

Then the function :func:`loop_package <pySecDec.loop_integral.loop_package>` is called. It will perform the algebraic sector decomposition steps and create a package containing the C++ code 
for the numerical evaluation. It will create a folder called box1L and allows to define parameters controlling the numerical part 
(for a complete list of possible options see  :func:`loop_package <pySecDec.loop_integral.loop_package>`).

.. code::

 loop_package(

 name = 'box1L',

 loop_integral = li,

 real_parameters = Mandelstam_symbols + mass_symbols,
 
 # complex_parameters are also possible

 # the highest order of the final epsilon expansion  
 requested_order = 0,

 # the optimization level to use in FORM (can be 0, 1, 2, 3)
 form_optimization_level = 2,

 # the WorkSpace parameter for FORM
 form_work_space = '100M',

 # the method to be used for the sector decomposition
 # valid values are ``iterative`` and ``geometric``
 decomposition_method = 'geometric',

 # whether or not to produce code to perform the contour deformation
 # if ``True``, it can still be deactivated later in the "config.hpp"
 # if ``False``, no code for the contour deformation is generated
 contour_deformation = True,

 )
 
Algebraic part and creation of the C++ library
---------------------------------------------- 

Running the python script  `box1L.py` 

.. code::

 $ python box1L.py
 
will create a folder with the name given in  `box1L.py`  ('box1L'),  which should contain the following files and subdirectories

.. code::

 box1L.hpp  integrate_box1L.cpp  box1L.pdf codegen  Makefile  Makefile.conf pylink README  src

in the folder 'box1L', typing 

.. code::

 $ make 
 
will create the libraries `libbox1L.a` and `box1L_pylink.so` which can be linked to an external program calling these integrals.
How to do this ``interactively`` or via a python script is explained in the next section. 
In ``standalone mode``, the C++ file `integrate_box1L.cpp` can be used to produce results for a certain kinematic point. In the latter, 
kinematic points can be specified by adapting the line
 
.. code::

     const std::vector<box1L::real_t> real_parameters = {9.,-0.1,0.3, 1.};
  
 
for the desired kinematics. In the above example, the values correspond to  `s=9,t=-0.1,s1=0.3, msq=1`, i.e. the same ordering is kept as in the lists Mandelstam_symbols = ['s','t','s1'],  mass_symbols = ['msq'] in the python input.

The commands 

.. code::

 $ make integrate_box1L
 $ ./integrate_box1L
 
will then evaluate the integral and print the result to the screen.



Interactive python interface
----------------------

There is also a python interface which allows for an interactive
evaluation of the integrals. 
We will use the 2-loop triangle example to explain how this works:

- first produce the code for the triangle by 

.. code::

 $ python triangle2L.py
 
- change to the directory triangle2L and type 

.. code::

 $ make 

- this produces, among other things,  the library  `triangle2L_pylink.so`. The latter can be called from within python. In order to do so,  ipython or python can be opened and the following commands can be entered interactively:

.. code::

>>> from __future__ import print_function
>>> from pySecDec.integral_interface import IntegralLibrary
>>> import sympy as sp
>>> # load c++ library
>>> triangle = IntegralLibrary('triangle2L_pylink.so')

- now the user can choose an integrator and define the settings for
  the numerical integration. A list of possible settings is given in :class:`pySecDec.integral_interface<pySecDec.integral_interface>`.

.. code::

>>> # choose integrator
>>> triangle.use_Vegas(flags=2,epsrel=1e-3,epsabs=1e-10) # ``flags=2`` means verbose --> see Cuba manual


- the numerical point at which the integral should be evaluated can be
  given as follows

.. code::

>>> # perform the integration for the numerical point s=0.9, msq=0.1
>>> str_integral_without_prefactor, str_prefactor,
>>> str_integral_with_prefactor = triangle(real_parameters=[.9,.1])

- the class *triangle* can take more parameters, for example

.. code::

>>> str_integral_with_prefactor = triangle(real_parameters=[.9,.1],number_of_presamples=1e+6,deformation_parameters_maximum = 0.5)
>>> #  (defaults: number_of_presamples = 100000, deformation_parameters_maximum = 1)

- further options for the contour deformation etc are listed under  :class:`pySecDec.integral_interface<pySecDec.integral_interface>`

- in addition, the output format can be specified:

.. code::

>>> # convert complex numbers from c++ to sympy notation
>>>  str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')
>>>  str_prefactor = str_prefactor.replace(',','+I*')
>>>  str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')

>>> # convert result to sympy expressions
>>>  integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
>>>  integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
>>>  prefactor = sp.sympify(str_prefactor)
>>>  integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
>>>  integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))

>>> # examples how to access individual orders
>>>  print('leading pole:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
>>>  print('subleading pole:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
>>>  print('finite part:', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (', integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')


- This will print the result in a format which is also easy to import into Mathematica. Examples for the above commands are also given in `integrate_triangle.py`. 

- How to loop over several kinematic points is shown in the example `multiple_points_example.py`.
