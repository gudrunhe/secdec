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

#define calIDerivatives "%(cal_I_derivatives)s"
#define functions "`calIDerivatives',%(functions)s"
CFunctions `functions';

#define decomposedPolynomialDerivatives "%(decomposed_polynomial_derivatives)s"
CFunctions `decomposedPolynomialDerivatives';

* Temporary functions and symbols for replacements in FORM
AutoDeclare CFunctions SecDecInternalfDUMMY;
AutoDeclare Symbols SecDecInternalsDUMMY;

* We generated logs in the subtraction and pack denominators
* and powers into a functions.
CFunctions log, SecDecInternalPow, SecDecInternalDenominator;

* We rewrite function calls as symbols
#Do function = {`functions',`decomposedPolynomialDerivatives',log,SecDecInternalPow,SecDecInternalDenominator}
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

* Specify and enumerate all occurring orders in python.
* Define the preprocessor variables
* `shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex''.
%(regulator_powers)s

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

* Define procedures to insert the dummy functions introduced in python and their derivatives.
#procedure insertCalI
  %(insert_cal_I_procedure)s
#endProcedure

#procedure insertOther
  %(insert_other_procedure)s
#endProcedure

#procedure insertDecomposed
  %(insert_decomposed_procedure)s
#endProcedure

* Define how deep functions to be inserted are nested.
#define insertionDepth "%(form_insertion_depth)i"

* Define the data type of the integrand container class (constructed in python).
#define integrandContainerType "%(sector_container_type)s"

* Define the initializer list for the integrand container class
* (constructed in python).
#define integrandContainerInitializer "%(sector_container_initializer)s"
