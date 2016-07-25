* The name of the loop integral
#define name "%(name)s"

* Whether or not we are producing code for contour deformation
#define contourDeformation "%(contour_deformation)i"

* number of integration variables
#define numIV "%(number_of_integration_variables)i"

* number of regulators
#define numReg "%(number_of_regulators)i"

#define integrationVariables "%(integration_variables)s"
#define realParameters "%(real_parameters)s"
#define complexParameters "%(complex_parameters)s"
#define regulators "%(regulators)s"
Symbols `integrationVariables'
        `realParameters'
        `complexParameters'
        `regulators';

* Define the imaginary unit in sympy notation.
Symbol I;

#define functions "%(functions)s"
CFunctions `functions';

*** TODO: disallow user variables with prefix "SecDecInternal" in python
* Temporary functions and symbols for replacements in FORM
AutoDeclare CFunctions SecDecInternalfDUMMY;
AutoDeclare Symbols SecDecInternalsDUMMY;

* TODO: How to determine which derivatives of the user input ``functions`` are needed? How to communicate it to the user? --> quick and dirty solution for a start:
AutoDeclare CFunctions d;
* TODO: remove the line above

* We generated logs in the subtraction and pack denominators into a function
CFunctions log, SecDecInternalDenominator;

* We rewrite function calls as symbols
#Do function = {`functions',log,SecDecInternalDenominator}
  AutoDeclare Symbols SecDecInternal`function'Call;
#EndDo

* We need labels for the code optimization
AutoDeclare Symbols SecDecInternalLabel;

* The integrand may be longer than FORM can read in one go.
* We use python to split the the expression if neccessary.
* Define a procedure that defines the "integrand" expression
#procedure defineExpansion
  Global expansion = SecDecInternalsDUMMYIntegrand;
  %(integrand_definition_procedure)s
#endProcedure

#define highestPoles "%(highest_regulator_poles)s"
#define numOrders "%(number_of_orders)s"
#define optimizationLevel "%(form_optimization_level)i"
#define stabilize "%(stabilize)i"

* Specify and enumerate all occurring orders in python.
* Define the preprocessor variables
* `shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex''.
%(regulator_powers)s

* The array of abbreviations
ExtraSymbols,array,SecDecInternalAbbreviation;

* Define two procedures to open and close a nested argument section
#procedure beginArgumentDepth(depth)
  #Do recursiveDepth = 1, `depth'
    Argument;
  #EndDo
#endProcedure
#procedure endArgumentDepth(depth)
  #Do recursiveDepth = 1, `depth'
    EndArgument;
  #EndDo
#endProcedure

* Define a procedure to insert the dummy functions introduced in python and their derivatives.
#procedure insert
  %(insert_procedure)s
#endProcedure

* Define how deep functions to be inserted are nested.
#define insertionDepth "%(form_insertion_depth)i"

* Define the data type of the integrand container class (constructed in python).
#define integrandContainerType "%(sector_container_type)s"

* Define the initializer list for the integrand container class
* (constructed in python).
#define integrandContainerInitializer "%(sector_container_initializer)s"
