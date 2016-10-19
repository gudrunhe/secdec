Getting started
===============

After installation, you should have a folder `examples` in your main pySecDec directory.
It contains examples for a one-loop box, `box1L.py`, a two-loop triangle, `triangle.py`, 
a two-loop box, `box2L.py`, and parametric function not related to loop integrals, `example_run_card_general.py`. 


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

 box1L.hpp  codegen  Makefile  Makefile.conf  README  src  test.cpp

in the folder 'box1L', typing 

.. code::

 $ make 
 
will create the libraries `libbox1L.a`,  `libbox1L.so` which can be linked to an external program calling these integrals.
In ``standalone mode``, the C++ file ``test.cpp`` can be used to produce results for a certain kinematic point. In the latter, 
kinematic points can be specified by adapting the line

.. TODO: or a list of kinematic points, read from a kinem.input file
 
.. code::

     const std::vector<box1L::real_t> real_parameters = {9.,-0.1,0.3, 1.};
  
 
for the desired kinematics. In the above example, the values correspond to  `s=9,t=-0.1,s1=0.3, msq=1`, i.e. the same ordering is kept as in the lists Mandelstam_symbols = ['s','t','s1'], 
 mass_symbols = ['msq'] in the python input.

The commands 

.. code::

 $ make test
 $ ./test
 
will then evaluate the integral and print the result to the screen.



