#-
Off statistics;

* define two procedures to open and close a nested argument section
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

#include sector`sectorID'.h
.global

Global expansion = `integrand';
.sort

* Enumerate the regulators
#$counter = 1;
#Do regulator = {`regulators'}
  #define regulator`$counter' "`regulator'"
  #$counter = $counter + 1;
#EndDo

* Enumerate the poles
#$counter = 1;
#Do highestPole = `highestPoles'
  #define highestPole`$counter' "`highestPole'"
  #$counter = $counter + 1;
#EndDo

* FORM is not good at handling negative powers (poles) --> multiply by the highest poles
#Do i = 1,`numReg'
  multiply `regulator`i''^`highestPole`i'';
#EndDo

* Bracket according to the regulators to separate the orders.
B `regulators';

* Optimize each order in epsilon separately.
* The orders to be processed are enumerated in python. The shifted power
* of each regulator is stored in the preprocessor variable
* `shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex''.
#Do shiftedOrderIndex = 1, `numOrders'
* clear previous step
  .store

* Calculate the (possibly negative) orders in the regulators
* They are used as multiindex for the function name.
  #Redefine cppOrder ""
  #Do regulatorIndex = 1, `numReg'
    #$absOfOrder = `shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex''-`highestPole`regulatorIndex'';

*   Since we are not allowed to have a "-" in c++ function names,
*   replace the "-" by an "n" if required
    #if `$absOfOrder' < 0
      #$absOfOrder = - $absOfOrder;
      #Redefine cppOrder "`cppOrder'n`$absOfOrder'"
    #else
      #Redefine cppOrder "`cppOrder'`$absOfOrder'"
    #endif

*   Separate the orders in the different regulators by underscores
    #if `regulatorIndex' != `numReg'
      #Redefine cppOrder "`cppOrder'_"
    #endif
  #EndDo

* We are writing a the c++ file "sector_`sectorID'_`cppOrder'.cpp"
* and the corresponding header "sector_`sectorID'_`cppOrder'.hpp".
* The header can already completely be written here:
  #write <sector_`sectorID'_`cppOrder'.hpp> "#ifndef __SecDec_include_guard_`name'_sector_`sectorID'_order_`cppOrder'"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#define __SecDec_include_guard_`name'_sector_`sectorID'_order_`cppOrder'"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#include <`name'/config.hpp>"
  #write <sector_`sectorID'_`cppOrder'.hpp> "namespace `name'"
  #write <sector_`sectorID'_`cppOrder'.hpp> "{"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  #if `name'_contour_deformation"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    IntegrandFunction"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  #else `name'_contour_deformation"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    DeformableIntegrandFunction"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  #endif"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  sector_`sectorID'__order_`cppOrder'__integrand;"
  #write <sector_`sectorID'_`cppOrder'.hpp> "}"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#endif"

* Open the namspace in which the sector is to be implemented
  #write <sector_`sectorID'_`cppOrder'.cpp> "#include <`name'/integrands/sector_`sectorID'_`cppOrder'.hpp>"
  #write <sector_`sectorID'_`cppOrder'.cpp> "namespace `name'"
  #write <sector_`sectorID'_`cppOrder'.cpp> "{"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  integrand_return_t sector_`sectorID'__order_`cppOrder'__integrand"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  ("
  #write <sector_`sectorID'_`cppOrder'.cpp> "    #if `name'_contour_deformation"
  #write <sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const integration_variables,"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    #else `name'_contour_deformation"
  #write <sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const integration_variables,"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    #endif"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    real_t const * const real_parameters,"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    complex_t const * const complex_parameters"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  )"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  {"

* extract the order in the regulators that we are about to process
  #$currentOrder = 1;
  #Do regulatorIndex = 1, `numReg'
    #$currentOrder = $currentOrder * `regulator`regulatorIndex''^`shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex'';
  #EndDo
  Local expression = expansion[$currentOrder];

* Explicitly insert the functions defined in python.
* {

* construct function argument for the left hand side of the `id` statement
* Example:
* Id f(x$sDUMMY1, y$sDUMMY2) = ...
  #redefine matchArg ""
  #$counter = 1;
  #Do var = {`integrationVariables', `regulators'}
    #If `$counter' != 1
      #redefine matchArg "`matchArg' , "
    #EndIf
    #redefine matchArg "`matchArg' `var'?$sDUMMY`$counter'"
    #$counter = $counter + 1;
  #EndDo

* perform innermost replacements first
  #Do minusDepth = - `insertionDepth' , 0
    #$depth = - `minusDepth';

*   Since we need intermediate ".sort" instructions, we cannot use the "repeat" environment.
*   The following construction is suggested in the FORM documentation.
    #Do i = 1,1

*     Cancel ratios of functions and wrap denominators into the function "SecDecInternalDenominator".
*     example: "U(x,y,z)/U(x,y,z)^2" --> "SecDecInternalDenominator(U(x,y,z))"
      #call beginArgumentDepth(`$depth')
        Denominators SecDecInternalDenominator;
        factarg SecDecInternalDenominator;
        chainout SecDecInternalDenominator;
        repeat Id fDUMMY?(?sDUMMY) * SecDecInternalDenominator(fDUMMY?(?sDUMMY)) = 1;
      #call endArgumentDepth(`$depth')
      .sort

*     expand logs
      #call beginArgumentDepth(`$depth')
        factarg log;
        repeat Id log(?head, sDUMMY1?, sDUMMY2?) = log(?head, sDUMMY1) + log(sDUMMY2);
        repeat Id log(sDUMMY1? ^ sDUMMY2?) = log(sDUMMY1) * sDUMMY2;
      #call endArgumentDepth(`$depth')
      .sort

*     some simplifications
      #call beginArgumentDepth(`$depth')
        Id log(1) = 0;
        repeat Id sDUMMY1? ^ sDUMMY2?neg_ = SecDecInternalDenominator(sDUMMY1) ^ (-sDUMMY2);
        repeat Id 1/sDUMMY? = SecDecInternalDenominator(sDUMMY);
        repeat Id sDUMMY? * SecDecInternalDenominator(sDUMMY?) = 1;
      #call endArgumentDepth(`$depth')
      .sort

      #Do functionForInsertion = {`functionsForInsertion'}

*       set dollar variables
        #call beginArgumentDepth(`$depth')
          if ( match(`functionForInsertion'(`matchArg')) ) redefine i "0";
        #call endArgumentDepth(`$depth')
        .sort

        #redefine replaceArg ""
        #$counter = 1;
        #Do var = {`integrationVariables', `regulators'}
          #If `$counter' != 1
            #redefine replaceArg "`replaceArg' , "
          #EndIf
          #redefine replaceArg "`replaceArg' `var',$sDUMMY`$counter'"
          #$counter = $counter + 1;
        #EndDo

        #redefine idArg ""
        #Do j = 1,`numIV'+`numReg'
          #If `j' != 1
            #redefine idArg "`idArg', "
          #EndIf
          #redefine idArg "`idArg'`$sDUMMY`j''"
        #EndDo

*       This "if" evaluates to true only if the expression above was matched.
*       If the expression above does not match, there is nothing to do.
        #If `i' == 0
          L replacement = (``functionForInsertion'')*replace_(`replaceArg');
          .sort
          drop replacement;
          #call beginArgumentDepth(`$depth')
            id `functionForInsertion'(`idArg') = replacement;
          #call endArgumentDepth(`$depth')
          .sort
        #EndIf

      #EndDo
    #EndDo
  #EndDo

* }

* Analytically cancel the subtraction terms to avoid numerical instabilities.
* We bring all terms that come with a "1/integration_variable" factor to a common
* denominator.
* {
  #If `stabilize'
    L denom = sDUMMYdenominator;
    #Do IV = {`integrationVariables'}
      #Do i = 1,1
        if ( match(SecDecInternalDenominator(`IV') * SecDecInternalDenominator(sDUMMY?!{`IV'}$arg)) ) redefine i "0";
        .sort
        #if `i' == 0
          multiply numerator($arg);
          Id numerator($arg) * sDUMMYdenominator = SecDecInternalDenominator($arg) * sDUMMYdenominator;
          Id numerator($arg) * SecDecInternalDenominator($arg) = 1;
          Id numerator($arg) = $arg;
          Id `FP' * SecDecInternalDenominator(`FP') = 1;
          .sort
        #EndIf
      #EndDo
    #EndDo
    Id sDUMMYdenominator = 1;
    .sort
    multiply denom;
    .sort
    drop denom;
    .sort
  #EndIf
* }

* Replace all function calls by symbols for simultaneous optimization.
* {
  Local toOptimize = sDUMMYtoOptimize;

  #Do function = {`functions',log,SecDecInternalDenominator}
    #$labelCounter = 0;

*   process innermost function calls first
    #Do minusDepth = - `insertionDepth' , 0
      #$depth = - `minusDepth';

*     Since we need intermediate ".sort" instructions, we cannot use the
*     "repeat" environment.
*     The following construction is suggested in the FORM documentation.

      #Do i = 1,1
*       set dollar variable
        #call beginArgumentDepth(`$depth')
          if ( match(`function'(?sDUMMY$args)) ) redefine i "0";
        #call endArgumentDepth(`$depth')
        .sort

*       The following "#if" evaluates to true only if there are logs or denominators left.
        #If `i' == 0

          #$labelCounter = $labelCounter + 1;

          L arguments = SecDecInternalLabel`function'Call`$labelCounter'Arg * fDUMMYarguments(`$args');

          #call beginArgumentDepth(`$depth')
            Id `function'(`$args') = `function'Call`$labelCounter';
          #call endArgumentDepth(`$depth')

          repeat Id fDUMMYarguments(sDUMMY?, ?otherArgs) = SecDecInternalLabel`function'Call`$labelCounter'Arg * (sDUMMY + fDUMMYarguments(?otherArgs));

*         Define `$argCounter' by loking at the term with the empty function "fDUMMYarguments"
          Id fDUMMYarguments * SecDecInternalLabel`function'Call`$labelCounter'Arg ^ sDUMMYexponent?$argCounter = 0;
          .sort
*         `$argCounter' as defined above is too big by one.
          #$argCounter = $argCounter - 1;

*         Add all arguments to top level polynomial for simultaneous optimization.
          Id sDUMMYtoOptimize = sDUMMYtoOptimize + arguments;

          #redefine numberOfArgs`function'Label`$labelCounter' "`$argCounter'"

          .sort
        #EndIf
      #EndDo
    #EndDo

    #redefine largestLabel`function' "`$labelCounter'"

  #EndDo

  Id sDUMMYtoOptimize = expression;
  .sort

  drop expression, arguments;
  .sort
* }

* Simultaneously optimize the integrand and all occuring function arguments.
  AB `integrationVariables', `realParameters', `complexParameters';
  Format O`optimizationLevel';
  .sort
  #optimize toOptimize

* Define the integration variables and parameters as c preprocessor variables
* (The integrand function in c takes them packed into an array).
* {
* "Format rational": Need the indices as integers.
  Format rational;

  #$counter = 0;
  #Do IV = {`integrationVariables'}
    #write <sector_`sectorID'_`cppOrder'.cpp> "#define `IV' integration_variables[`$counter']"
    #$counter = $counter + 1;
  #EndDo
  #$counter = 0;
  #Do RP = {`realParameters'}
    #write <sector_`sectorID'_`cppOrder'.cpp> "#define `RP' real_parameters[`$counter']"
    #$counter = $counter + 1;
  #EndDo
  #$counter = 0;
  #Do CP = {`complexParameters'}
    #write <sector_`sectorID'_`cppOrder'.cpp> "#define `CP' complex_parameters[`$counter']"
    #$counter = $counter + 1;
  #EndDo
* }

* Processing denominators in FORM is easiest if packed into a function.
* Define that function as c preprocessor macro.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalDenominator(x) 1./(x)"

* Define "SecDecInternalAbbreviation[0]" as c preprocessor variable "result".
* Since FORM does not use "SecDecInternalAbbreviation[0]", we can use it.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define result SecDecInternalAbbreviation[0]"

* write Abbreviations in c format
  Format float 20;
  Format C;
  #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t"
  #write <sector_`sectorID'_`cppOrder'.cpp> "SecDecInternalAbbreviation[`optimmaxvar_' + 1];"
  #write <sector_`sectorID'_`cppOrder'.cpp> "%%O"

* Replace all function calls by symbols for simultaneous optimization.
* {
  #Do function = {`functions',log,SecDecInternalDenominator}
    #Do callIndex = 1, `largestLabel`function''
      B label`function'Call`callIndex'Arg;
      .sort
      #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
        L arg`argIndex' = toOptimize[label`function'Call`callIndex'Arg ^ `numberOfArgs`function'Label`callIndex'];
      #EndDo
      .sort
      #write <sector_`sectorID'_`cppOrder'.cpp> ""
      #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t `function'Call`callIndex' ="
      #write <sector_`sectorID'_`cppOrder'.cpp> "`function'("
      #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
        #If `argIndex' == `numberOfArgs`function'Label`callIndex''
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%E"  arg`argIndex'
        #Else
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%E," arg`argIndex'
        #EndIf
        drop arg`argIndex';
      #EndDo
      #write <sector_`sectorID'_`cppOrder'.cpp> ");" currentExpr
      #write <sector_`sectorID'_`cppOrder'.cpp> ""
      multiply replace_(label`function'Call`callIndex'Arg, 0);
      .sort
    #EndDo
  #EndDo
* }

* write the integrand
  #write <sector_`sectorID'_`cppOrder'.cpp> ""
  #write <sector_`sectorID'_`cppOrder'.cpp> "result = %%e" toOptimize(result)
  #write <sector_`sectorID'_`cppOrder'.cpp> "return result;"
  #write <sector_`sectorID'_`cppOrder'.cpp> ""

* undefine the c preprocessor macros
  #Do IV = {`integrationVariables'}
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef `IV'"
  #EndDo
  #Do RP = {`realParameters'}
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef `RP"
  #EndDo
  #Do CP = {`complexParameters'}
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef `CP'"
  #EndDo
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalDenominator"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef result"

* Close the c++ function and namespaces
  #write <sector_`sectorID'_`cppOrder'.cpp> "  };"
  #write <sector_`sectorID'_`cppOrder'.cpp> "};"
#EndDo

******* TODO: continue adaptation of this file below

******* TODO: how to initialize the multivariate series ---> finish python's algebraic part before bothering with the syntax here
* Write a c++ header that collects all the functions in a Series
* "Format rational": Need the indices as integers.
Format rational;
#write <sector_`sectorID'.hpp> "#ifndef __SecDec_include_guard_`name'_sector_`sectorID'"
#write <sector_`sectorID'.hpp> "#define __SecDec_include_guard_`name'_sector_`sectorID'"
#write <sector_`sectorID'.hpp> "#include <`name'/util/series.hpp>"
#Do shiftedOrder = 0, `numOrders'
* Calculate the (possibly negative) order in the regulator
  #$order = `shiftedOrder'-`highestPole';
* Since we are not allowed to have a "-" in c++ function names,
* replace the "-" by an "n" if required
  #if `$order' < 0
    #$absOfOrder = - $order;
    #Redefine cppOrder "n`$absOfOrder'"
  #else
    #$absOfOrder = + $order;
    #Redefine cppOrder "`$order'"
  #endif
* include correspoinding c++ header file
#write <sector_`sectorID'.hpp> "#include <`name'/integrands/sector_`sectorID'_`cppOrder'.hpp>"
#EndDo
#write <sector_`sectorID'.hpp> "#if `name'_contour_deformation"
#write <sector_`sectorID'.hpp> "#include <`name'/integrands/contour_deformation_sector_`sectorID'.hpp>"
#write <sector_`sectorID'.hpp> "#endif"
#write <sector_`sectorID'.hpp> "namespace `name'"
#write <sector_`sectorID'.hpp> "{"
#write <sector_`sectorID'.hpp> "  Series<integrand_t> integrand_of_sector_`sectorID'"
#write <sector_`sectorID'.hpp> "  ("
#write <sector_`sectorID'.hpp> "      /*order_min*/ -`highestPole',"
#write <sector_`sectorID'.hpp> "      /*truncated_below*/ false,"
#write <sector_`sectorID'.hpp> "      /*order_max*/ `numOrders' - (`highestPole'),"
#write <sector_`sectorID'.hpp> "      /*truncated_above*/ true,"
#write <sector_`sectorID'.hpp> "      {"
#Do shiftedOrder = 0, `numOrders'
* Calculate the (possibly negative) order in the regulator
  #$order = `shiftedOrder'-`highestPole';
* Since we are not allowed to have a "-" in c++ function names,
* replace the "-" by an "n" if required
  #if `$order' < 0
    #$absOfOrder = - $order;
    #Redefine cppOrder "n`$absOfOrder'"
  #else
    #$absOfOrder = + $order;
    #Redefine cppOrder "`$order'"
  #endif
* add corresponding function to the Series
#write <sector_`sectorID'.hpp> "{"
#write <sector_`sectorID'.hpp> "/*number of Feynman parameters*/ `NumFP',"
#write <sector_`sectorID'.hpp> "sector_`sectorID'__order_`cppOrder'__integrand,"
#write <sector_`sectorID'.hpp> "#if `name'_contour_deformation"
#write <sector_`sectorID'.hpp> "sector_`sectorID'__contourdef,"
#write <sector_`sectorID'.hpp> "sector_`sectorID'__F"
#write <sector_`sectorID'.hpp> "#endif"
#write <sector_`sectorID'.hpp> "},"
#EndDo
#write <sector_`sectorID'.hpp> "      }"
#write <sector_`sectorID'.hpp> "  );"
#write <sector_`sectorID'.hpp> "};"
#write <sector_`sectorID'.hpp> "#endif"

.end
