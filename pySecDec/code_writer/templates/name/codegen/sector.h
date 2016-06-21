* The name of the loop integral
#define name "%(name)s"

* number of integration variables
#define numIV "%(number_of_integration_variables)i"
#define numReg "%(number_of_regulators)i"

#define integrationVariables "%(integration_variables)s"
#define realParameters "%(real_parameters)s"
#define complexParameters "%(complex_parameters)s"
#define regulators "%(regulators)s"
Symbols `integrationVariables'
        `realParameters'
        `complexParameters'
        `regulators';

#define functions "%(functions)s"
CFunctions `functions';

*** TODO: add prefix "SecDecInternal" to all internal variables and disallow user variables with that prefix
* Temporary functions and symbols for replacements in FORM
AutoDeclare CFunctions fDUMMY;
AutoDeclare Symbols sDUMMY;

* TODO: How to determine which derivatives of the user input ``functions`` are needed? How to communicate it to the user? --> quick and dirty solution for a start:
AutoDeclare CFunctions d;
* TODO: remove the line above

* We generated logs in the subtraction and pack denominators into a function
CFunctions log, SecDecInternalDenominator;

* We rewrite function calls as symbols
AutoDeclare Symbols SecDecInternalCall, SecDecInternalDenominatorCall;
#Do function = {`functions',log}
  AutoDeclare Symbols `function'Call;
#EndDo

* We need labels for the code optimization
AutoDeclare Symbols SecDecInternalLabel;

#define integrand "%(integrand)s"

#define highestPoles "%(highest_regulator_poles)s"
#define numOrders "%(number_of_orders)s"
#define optimizationLevel "%(form_optimization_level)i"
#define stabilize "%(stabilize)i"

* Specify and enumerate all occurring orders in python.
* Define the preprocessor variables
* `shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex''.
%(regulator_powers)s

* The array of abbreviations
ExtraSymbols,array,abbreviation;

* Define the dummy functions introduced in python and their derivatives.
#define functionsForInsertion "%(functions_for_insertion)s"
%(function_definitions)s

* Define how deep functions to be inserted are nested.
#define insertionDepth "%(form_insertion_depth)i"
