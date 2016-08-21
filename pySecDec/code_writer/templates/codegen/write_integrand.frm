#-
Off statistics;

* Define two general procedures that write c++ code to define and undefine
* c++ preprocessor varibables accessing a c++ array.
#procedure cppDefine(?FORMNames,cppArrayName)
  #$counter = 0;
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <sector_`sectorID'_`cppOrder'.cpp> "#define `varname' `cppArrayName'[`$counter']#@SecDecInternalNewline@#"
      #$counter = $counter + 1;
    #EndIf
  #EndDo
#endProcedure
#procedure cppUndefine(?FORMNames)
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <sector_`sectorID'_`cppOrder'.cpp> "#undef `varname'#@SecDecInternalNewline@#"
    #EndIf
  #EndDo
#endProcedure

* Define the same two general procedures for the contour deformation
* c++ preprocessor varibables accessing a c++ array.
#procedure cppDefineContourdef(?FORMNames,cppArrayName)
  #$counter = 0;
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define `varname' `cppArrayName'[`$counter']#@SecDecInternalNewline@#"
      #$counter = $counter + 1;
    #EndIf
  #EndDo
#endProcedure
#procedure cppUndefineContourdef(?FORMNames)
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#undef `varname'#@SecDecInternalNewline@#"
    #EndIf
  #EndDo
#endProcedure

#include sector`sectorID'.h
#If `contourDeformation'
  #include contour_deformation_sector`sectorID'.h
#EndIf
.global

#call defineExpansion
.sort

* Enumerate the regulators
#$counter = 1;
#Do regulator = {`regulators'}
  #define regulator`$counter' "`regulator'"
  #$counter = $counter + 1;
#EndDo

* Enumerate the poles
#$counter = 1;
#Do highestPole = {`highestPoles',}
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
  #write <sector_`sectorID'_`cppOrder'.hpp> "#ifndef `name'_codegen_sector_`sectorID'_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#define `name'_codegen_sector_`sectorID'_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#include \"`name'.hpp\"#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "{#@SecDecInternalNewline@#"
  #If `contourDeformation'
    #write <sector_`sectorID'_`cppOrder'.hpp> "    secdecutil::SectorContainerWithDeformation<real_t, complex_t>::DeformableIntegrandFunction#@SecDecInternalNewline@#"
  #Else
    #write <sector_`sectorID'_`cppOrder'.hpp> "    secdecutil::SectorContainerWithoutDeformation<real_t, complex_t, integrand_return_t>::IntegrandFunction#@SecDecInternalNewline@#"
  #EndIf
  #write <sector_`sectorID'_`cppOrder'.hpp> "  sector_`sectorID'_order_`cppOrder'_integrand;#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "}#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#endif#@SecDecInternalNewline@#"

* Open the namspace in which the sector is to be implemented
  #write <sector_`sectorID'_`cppOrder'.cpp> "#include \"sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "{#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  integrand_return_t sector_`sectorID'_order_`cppOrder'_integrand#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  (#@SecDecInternalNewline@#"
  #If `contourDeformation'
    #write <sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const integration_variables,#@SecDecInternalNewline@#"
  #Else
    #write <sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const integration_variables,#@SecDecInternalNewline@#"
  #EndIf
  #write <sector_`sectorID'_`cppOrder'.cpp> "    real_t const * const real_parameters,#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    complex_t const * const complex_parameters#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  )#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "  {#@SecDecInternalNewline@#"

* extract the order in the regulators that we are about to process
  #$currentOrder = 1;
  #Do regulatorIndex = 1, `numReg'
    #$currentOrder = $currentOrder * `regulator`regulatorIndex''^`shiftedRegulator`regulatorIndex'PowerOrder`shiftedOrderIndex'';
  #EndDo
  Local expression = expansion[$currentOrder];

* Expand logs *BEFORE* insertions and only at top level in order to avoid
* introducing log(<negative real>).
  factarg log;
  repeat Id log(?head, SecDecInternalsDUMMY1?, SecDecInternalsDUMMY2?) = log(?head, SecDecInternalsDUMMY1) + log(SecDecInternalsDUMMY2);
  repeat Id log(SecDecInternalsDUMMY1? ^ SecDecInternalsDUMMY2?) = log(SecDecInternalsDUMMY1) * SecDecInternalsDUMMY2;
  .sort

* Explicitly insert the functions defined in python.
* {

  #Do depth = 0, `insertionDepth'
*   Cancel ratios of functions and wrap denominators into the function "SecDecInternalDenominator".
*   example: "U(x,y,z)/U(x,y,z)^2" --> "SecDecInternalDenominator(U(x,y,z))"
    #call beginArgumentDepth(`depth')
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;
    #call endArgumentDepth(`depth')
    .sort

    #call beginArgumentDepth(`depth')
      #call insert
    #call endArgumentDepth(`depth')
    .sort

*   some simplifications
    #call beginArgumentDepth(`depth')
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      Id log(1) = 0;
      repeat Id SecDecInternalsDUMMY1? ^ SecDecInternalsDUMMY2?neg_ = SecDecInternalDenominator(SecDecInternalsDUMMY1) ^ (-SecDecInternalsDUMMY2);
      repeat Id 1/SecDecInternalsDUMMY? = SecDecInternalDenominator(SecDecInternalsDUMMY);
      repeat Id SecDecInternalsDUMMY? * SecDecInternalDenominator(SecDecInternalsDUMMY?) = 1;
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;
      repeat Id SecDecInternalDenominator(SecDecInternalsDUMMY?number_) = 1/SecDecInternalsDUMMY;
    #call endArgumentDepth(`depth')
    .sort
  #EndDo

* }

* simplify again
  #Do depth = 0, `insertionDepth'
    #call beginArgumentDepth(`depth')
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      Id log(1) = 0;
      repeat Id SecDecInternalsDUMMY1? ^ SecDecInternalsDUMMY2?neg_ = SecDecInternalDenominator(SecDecInternalsDUMMY1) ^ (-SecDecInternalsDUMMY2);
      repeat Id 1/SecDecInternalsDUMMY? = SecDecInternalDenominator(SecDecInternalsDUMMY);
      repeat Id SecDecInternalsDUMMY? * SecDecInternalDenominator(SecDecInternalsDUMMY?) = 1;
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;
      repeat Id SecDecInternalDenominator(SecDecInternalsDUMMY?number_) = 1/SecDecInternalsDUMMY;
    #call endArgumentDepth(`depth')
    .sort
  #EndDo

* Analytically cancel the subtraction terms to avoid numerical instabilities.
* We bring all terms that come with a "1/integration_variable" factor to a common
* denominator.
* {
  #If `stabilize'
    L denom = SecDecInternalsDUMMYdenominator;
    #Do IV = {`integrationVariables'}
      #Do i = 1,1
        if ( match(SecDecInternalDenominator(`IV') * SecDecInternalDenominator(SecDecInternalsDUMMY?!{`IV'}$arg)) ) redefine i "0";
        .sort
        #if `i' == 0
          multiply SecDecInternalfDUMMYNumerator($arg);
          Id SecDecInternalfDUMMYNumerator($arg) * SecDecInternalsDUMMYdenominator = SecDecInternalDenominator($arg) * SecDecInternalsDUMMYdenominator;
          Id SecDecInternalfDUMMYNumerator($arg) * SecDecInternalDenominator($arg) = 1;
          Id SecDecInternalfDUMMYNumerator($arg) = $arg;
          Id `IV' * SecDecInternalDenominator(`IV') = 1;
          .sort
        #EndIf
      #EndDo
    #EndDo
    Id SecDecInternalsDUMMYdenominator = 1;
    .sort
    multiply denom;
    .sort
    drop denom;
    .sort
  #EndIf
* }

* Replace all function calls by symbols for simultaneous optimization.
* {
  Local toOptimize = SecDecInternalsDUMMYtoOptimize;

  #Do function = {`functions',log,SecDecInternalDenominator}
    #$labelCounter = 0;

    #Do depth = 0, `insertionDepth'

*     Since we need intermediate ".sort" instructions, we cannot use the
*     "repeat" environment.
*     The following construction is suggested in the FORM documentation.

      #Do i = 1,1
*       set dollar variable
        #call beginArgumentDepth(`depth')
          if ( match(`function'(?SecDecInternalsDUMMY$args)) ) redefine i "0";
        #call endArgumentDepth(`depth')
        .sort

*       The following "#if" evaluates to true only if there are logs or denominators left.
        #If `i' == 0

          #$labelCounter = $labelCounter + 1;

          L arguments = SecDecInternalfDUMMYarguments(`$args');

          #call beginArgumentDepth(`depth')
            Id `function'(`$args') = SecDecInternal`function'Call`$labelCounter';
          #call endArgumentDepth(`depth')

          repeat Id SecDecInternalfDUMMYarguments(SecDecInternalsDUMMY?, ?otherArgs) = SecDecInternalLabel`function'Call`$labelCounter'Arg * (SecDecInternalsDUMMY + SecDecInternalfDUMMYarguments(?otherArgs));

*         Define `$argCounter' by loking at the term with the empty function "SecDecInternalfDUMMYarguments"
          Id SecDecInternalfDUMMYarguments * SecDecInternalLabel`function'Call`$labelCounter'Arg ^ SecDecInternalsDUMMYexponent?$argCounter = 0;
          .sort

*         Add all arguments to top level polynomial for simultaneous optimization.
          Id SecDecInternalsDUMMYtoOptimize = SecDecInternalsDUMMYtoOptimize + arguments;

          #redefine numberOfArgs`function'Label`$labelCounter' "`$argCounter'"

          .sort
        #EndIf
      #EndDo
    #EndDo

    #redefine largestLabel`function' "`$labelCounter'"

  #EndDo

  Id SecDecInternalsDUMMYtoOptimize = expression;
  .sort

  drop expression, arguments;
  .sort
* }

* translate sympy's imaginary unit to FORM's imaginary unit
multiply replace_(I,i_);
.sort

* Find and count the occurring integration variables.
* {
  #redefine occurringIntegrationVariables ""
  #redefine absentIntegrationVariables ""
  #$counterOccur = 0;
  #$counterAbsent = 0;
  #$currentIVIndex = -1;

  #Do IV = {`integrationVariables'}

    #$currentIVIndex = $currentIVIndex + 1;

    #redefine IVOccurs "0"
    if ( occurs(`IV') ) redefine IVOccurs "1";
    .sort

    #If `IVOccurs'
      #$counterOccur = $counterOccur + 1;
      #If `$counterOccur' == 1
        #redefine occurringIntegrationVariables "`IV'"
        #redefine occurringIntegrationVariableIndices "`$currentIVIndex'"
      #Else
        #redefine occurringIntegrationVariables "`occurringIntegrationVariables',`IV'"
        #redefine occurringIntegrationVariableIndices "`occurringIntegrationVariableIndices',`$currentIVIndex'"
      #EndIf
    #Else
      #$counterAbsent = $counterAbsent + 1;
      #If `$counterAbsent' == 1
        #redefine absentIntegrationVariables "`IV'"
        #redefine absentIntegrationVariableIndices "`$currentIVIndex'"
      #Else
        #redefine absentIntegrationVariables "`absentIntegrationVariables',`IV'"
        #redefine absentIntegrationVariableIndices "`absentIntegrationVariableIndices',`$currentIVIndex'"
      #EndIf
    #EndIf

  #EndDo

  #redefine numOccurringIVOrder`shiftedOrderIndex' "`$counterOccur'"
* }

* Specify the occurring deformation parameters.
  #If `contourDeformation'
    #redefine occurringDeformationParameters ""
    #$counter = 1;
    #Do var = {`occurringIntegrationVariableIndices',}
      #If x`var' != x
        #If `$counter' != 1
           #redefine occurringDeformationParameters "`occurringDeformationParameters',"
        #Endif
        #redefine occurringDeformationParameters "`occurringDeformationParameters'SecDecInternalLambda`var'"
        #$counter = $counter + 1;
      #EndIf
    #EndDo
  #EndIf


* Simultaneously optimize the integrand and all occurring function arguments.
  AB `integrationVariables', `realParameters', `complexParameters';
  Format O`optimizationLevel';
  .sort
  #optimize toOptimize

* Define the integration variables and parameters as c preprocessor variables
* (The integrand function in c takes them packed into an array).
* {
* "Format rational": Need the indices as integers.
  Format rational;

* call the general procedure to write the corresponding c++ code define in the beginning of this file
  #call cppDefine(`occurringIntegrationVariables',integration_variables)
  #call cppDefine(`realParameters',real_parameters)
  #call cppDefine(`complexParameters',complex_parameters)
* }

* Processing denominators in FORM is easiest if packed into a function.
* Define that function as c preprocessor macro.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalDenominator(x) 1./(x)#@SecDecInternalNewline@#"

* Define "SecDecInternalAbbreviation[0]" as c preprocessor variable "result".
* Since FORM does not use "SecDecInternalAbbreviation[0]", we can use it.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define result SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* write Abbreviations in c format
  Format float 20;
  Format C;
  #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* Replace all function calls by symbols for simultaneous optimization.
* {
  #Do function = {`functions',log,SecDecInternalDenominator}
    #Do callIndex = 1, `largestLabel`function''
      B SecDecInternalLabel`function'Call`callIndex'Arg;
      .sort
      #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
        L arg`argIndex' = toOptimize[SecDecInternalLabel`function'Call`callIndex'Arg ^ `argIndex'];
      #EndDo
      .sort
      #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternal`function'Call`callIndex' = "
      #write <sector_`sectorID'_`cppOrder'.cpp> "`function'("
      #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
        #If `argIndex' == `numberOfArgs`function'Label`callIndex''
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%E"  arg`argIndex'
        #Else
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%E," arg`argIndex'
        #EndIf
        drop arg`argIndex';
      #EndDo
      #write <sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#" currentExpr
      multiply replace_(SecDecInternalLabel`function'Call`callIndex'Arg, 0);
      .sort
    #EndDo
  #EndDo
* }

* write the integrand
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "result = %%e#@SecDecInternalNewline@#" toOptimize(result)
  #write <sector_`sectorID'_`cppOrder'.cpp> "return result;#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables')
  #call cppUndefine(`realParameters')
  #call cppUndefine(`complexParameters')
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalDenominator#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef result#@SecDecInternalNewline@#"

* Close the c++ function and namespaces
  #write <sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* write the contour deformation if required

* Must define preprocessor variable `F' even if we do not need contour deformation
* since FORM crashes otherwise.
  #If `contourDeformation' == 0
    #redefine F ""
  #EndIf

  #If `contourDeformation'
      #clearoptimize
      drop;
      Format normal;
      Format rational;
      .sort

*     We are writing the c++ file "contour_deformation_sector_`sectorID'_`cppOrder'.cpp"
*     and the corresponding header "contour_deformation_sector_`sectorID'_`cppOrder'.hpp".
*     The header can already completely be written here:
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#ifndef `name'_codegen_contour_deformation_sector_`sectorID'_order_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#define `name'_codegen_contour_deformation_sector_`sectorID'_order_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include \"`name'.hpp\"#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_complex_math.h>#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_linalg.h>#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "namespace `name'#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "{#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "  secdecutil::SectorContainerWithDeformation<real_t, complex_t>::ContourDeformationFunction sector_`sectorID'_order_`cppOrder'_contour_deformation;#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "  secdecutil::SectorContainerWithDeformation<real_t, complex_t>::DeformableIntegrandFunction sector_`sectorID'_order_`cppOrder'_contour_deformation_polynomial;#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "};#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#endif#@SecDecInternalNewline@#"

*     Open the namespace and the function in the corresponding .cpp file
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "namespace `name'#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "{#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  secdecutil::integral_transformation_t<complex_t> sector_`sectorID'_order_`cppOrder'_contour_deformation#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  (#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const integration_variables,#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const complex_parameters,#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const deformation_parameters#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  )#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  {#@SecDecInternalNewline@#"

*     The following procedure creates a local expression
*     with the integral transformation (the "contour deformation")
*     and its Jacobian matrix suitable for simultaneous optimization.
      Local contourdef = SecDecInternalsDUMMYContourdefExpression + SecDecInternalsDUMMYContourdefAppendix;
      #call insertContourdefExpression
      .sort

*     Implement commuting derivatives for the contour deformation polynomial F;
*     i.e. replace ddFdjdi --> ddFdidj for i<=j.
*     This is neccessary because pySecDec's algebra module does not
*     imply commuting derivatives for general functions.
      #Do i = 0,`$numIVMinusOne'
        #Do j = `i',`$numIVMinusOne'
          Id dd`F'd`j'd`i'(?args) = dd`F'd`i'd`j'(?args);
        #EndDo
      #EndDo
      .sort

*     We must take only the real part of `F' and its derivatives.
*     Collect all these and replace calls by symbols.
      Id `F'(?args) = SecDecInternal`F'Call;
      Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                              `F'(`integrationVariables',`regulators')*SecDecInternalLabel`F';
      #Do idx1 = {`occurringIntegrationVariableIndices'}
        Id d`F'd`idx1'(?args) = SecDecInternald`F'd`idx1'Call;
        Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                                d`F'd`idx1'(`integrationVariables',`regulators')*SecDecInternalLabeld`F'd`idx1';
        #Do idx2 = {`occurringIntegrationVariableIndices'}
          #If `idx1' <= `idx2'
            Id dd`F'd`idx1'd`idx2'(?args) = SecDecInternaldd`F'd`idx1'd`idx2'Call;
            Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                                    dd`F'd`idx1'd`idx2'(`integrationVariables',`regulators')*SecDecInternalLabeldd`F'd`idx1'd`idx2';
          #EndIf
        #EndDo
      #EndDo
      Id SecDecInternalsDUMMYContourdefAppendix = 0;

*     Remove derivatives of `F' by absent integration variables
      #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
        #Do idx1 = {`absentIntegrationVariableIndices',}
          #If x`idx1' != x
            Id d`F'd`idx1'(?args) = 0;
            #Do idx2 = {`absentIntegrationVariableIndices',}
              #If x`idx2' != x
                Id dd`F'd`idx1'd`idx2'(?args) = 0;
              #EndIf
            #EndDo
          #EndIf
        #EndDo
      #EndIf

*     define the argument for "replace_" that sets all regulators to zero
      #redefine nullifyRegulators ""
      #$counter = 1;
      #Do var = {`regulators'}
        #If `$counter' != 1
          #redefine nullifyRegulators "`nullifyRegulators' , "
        #EndIf
        #redefine nullifyRegulators "`nullifyRegulators' `var',0"
        #$counter = $counter + 1;
      #EndDo

*     set all regulators and "SecDecInternalsDUMMYToOptimize" to zero
      multiply replace_(`nullifyRegulators' , SecDecInternalsDUMMYToOptimize,0);
      .sort

      Format rational;
      Format Normal;

*     define the argument for "replace_" that sets all absent integration variables to zero
      #redefine nullifyAbsentIVs ""
      #$counter = 1;
      #Do var = {`absentIntegrationVariables',}
        #If x`var' != x
          #If `$counter' != 1
            #redefine nullifyAbsentIVs "`nullifyAbsentIVs' , "
          #EndIf
          #redefine nullifyAbsentIVs "`nullifyAbsentIVs' `var',0"
          #$counter = $counter + 1;
        #EndIf
      #EndDo

*     Remove absent integration variables
      #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
        multiply replace_(`nullifyAbsentIVs');
        .sort
      #EndIf

*     Explicitly insert `F' and its derivatives
      #call insert

*     translate sympy's imaginary unit to FORM's imaginary unit
      multiply replace_(I,i_);
      .sort

*     Define the integration variables, the real parameters, the complex parameters,
*     and the deformation parameters as c preprocessor variables.
*     (The c function takes them packed into an array).
*     "Format rational": Need the indices as integers.
      Format rational;
      #call cppDefineContourdef(`occurringIntegrationVariables',integration_variables)
      #call cppDefineContourdef(`realParameters',real_parameters)
      #call cppDefineContourdef(`complexParameters',complex_parameters)
      #call cppDefineContourdef(`occurringDeformationParameters',deformation_parameters)

*     optimize
      AB `integrationVariables', `realParameters', `complexParameters';
      Format O`optimizationLevel';
      .sort
      #optimize contourdef

*     Since FORM does not use "abbreviation[0]", we can use it as temporary variable.
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

*     write the optimization symbols
      Format float 20;
      Format C;
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

*     define the symbols "SecDecInternal...Call" as c++ variables
*     {

      Format rational;
      Format C;
      Bracket SecDecInternalLabel`F';
      .sort
      L expr = contourdef[SecDecInternalLabel`F'];
      .sort
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(tmp)
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternal`F'Call = tmp.real();#@SecDecInternalNewline@#"

      #Do idx1 = {`occurringIntegrationVariableIndices'}
        Bracket SecDecInternalLabeld`F'd`idx1';
        .sort
        L expr = contourdef[SecDecInternalLabeld`F'd`idx1'];
        .sort
        #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(tmp)
        #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternald`F'd`idx1'Call = tmp.real();#@SecDecInternalNewline@#"

        #Do idx2 = {`occurringIntegrationVariableIndices'}
          #If `idx1' <= `idx2'
        Bracket SecDecInternalLabeldd`F'd`idx1'd`idx2';
        .sort
        L expr = contourdef[SecDecInternalLabeldd`F'd`idx1'd`idx2'];
        .sort
        #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(tmp)
        #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternaldd`F'd`idx1'd`idx2'Call = tmp.real();#@SecDecInternalNewline@#"
          #EndIf
        #endDo
      #endDo

*     }

*     write the transformation
      Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
      .sort
      Format float 20;
      Format C;
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "std::vector<complex_t> transformed_parameters(`numIV');#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex *Jacobian = #@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    gsl_matrix_complex_alloc(`numOccurringIVOrder`shiftedOrderIndex'',`numOccurringIVOrder`shiftedOrderIndex'');#@SecDecInternalNewline@#"
      #$i = -1;
      #Do idx1 = {`occurringIntegrationVariableIndices'}
        #$idx1plus1DollarVar = `idx1' + 1;
        Format rational;
        #redefine idx1plus1 "`$idx1plus1DollarVar'"
        Format float 20;
        Format C;
        #$i = $i + 1;
        Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
        .sort
        L expr = contourdef[SecDecInternalLabelTransformation^`idx1plus1'];
        .sort
        #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "transformed_parameters[`$i'] = %%e" expr(transformed_parameters[`$i'])
        #$j = -1;
        #Do idx2 = {`occurringIntegrationVariableIndices'}
          #$idx2plus1DollarVar = `idx2' + 1;
          Format rational;
          #redefine idx2plus1 "`$idx2plus1DollarVar'"
          Format float 20;
          Format C;
          #$j = $j + 1;
          Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
          .sort
          L expr = contourdef[SecDecInternalLabelJacobianMatrixI^`idx1plus1' * SecDecInternalLabelJacobianMatrixJ^`idx2plus1'];
          .sort
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e" expr(tmp)
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex_set#@SecDecInternalNewline@#"
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "(#@SecDecInternalNewline@#"
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    Jacobian,#@SecDecInternalNewline@#"
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    `$i', `$j',#@SecDecInternalNewline@#"
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    {tmp.real(), tmp.imag()}#@SecDecInternalNewline@#"
          #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
        #EndDo
      #EndDo
*     calculate the determinant numerically using the gsl
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_permutation * permutation = gsl_permutation_alloc (`numOccurringIVOrder`shiftedOrderIndex'');#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "int signum;#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_linalg_complex_LU_decomp (Jacobian, permutation, &signum);#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_complex det = gsl_linalg_complex_LU_det(Jacobian, signum);#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = {GSL_REAL(det), GSL_IMAG(det)};#@SecDecInternalNewline@#"

*     free manually allocated memory
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_permutation_free(permutation);#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex_free(Jacobian);#@SecDecInternalNewline@#"

      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "return secdecutil::integral_transformation_t<complex_t>#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    {#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "        std::move(transformed_parameters),#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "        std::move(tmp)#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    };#@SecDecInternalNewline@#"

*     Close the function in the c++ file
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"

*     Optmimize the 2nd Symnanzik polynomial "F"
      #clearoptimize
      drop contourdef, expr;
      L expressionF = `F'(`integrationVariables',`regulators')*replace_(`nullifyRegulators');
      .sort
      #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
        multiply replace_(`nullifyAbsentIVs');
        .sort
      #EndIf
      #call insert
      .sort
      #optimize expressionF

      Format float 20;
      Format C;

*     Write the function to optimize the contour deformation parameters
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  integrand_return_t sector_`sectorID'_order_`cppOrder'_contour_deformation_polynomial#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  (#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const integration_variables,#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const complex_parameters#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  )#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  {#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e" expressionF(tmp)
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "return tmp;#@SecDecInternalNewline@#"

*     undefine the c preprocessor macros
      #call cppUndefineContourdef(`occurringIntegrationVariables')
      #call cppUndefineContourdef(`realParameters')
      #call cppUndefineContourdef(`complexParameters')
      #call cppUndefineContourdef(`occurringDeformationParameters')
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"

*     Close the function and the namespace in the c++ file
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

  #EndIf

#EndDo


* Here, all integrand functions are written to the hard disc.
* We still need a header file that collects the whole sector.

* clear last step
.store

* "Format rational": Need the indices as integers.
Format rational;

* open the c++ include guard
#write <sector_`sectorID'.hpp> "#ifndef `name'_codegen_sector_`sectorID'_hpp_included#@SecDecInternalNewline@#"
#write <sector_`sectorID'.hpp> "#define `name'_codegen_sector_`sectorID'_hpp_included#@SecDecInternalNewline@#"

* include series class
#write <sector_`sectorID'.hpp> "#include <secdecutil/series.hpp>#@SecDecInternalNewline@#"

#write <sector_`sectorID'.hpp> "#@SecDecInternalNewline@#"

#Do shiftedOrderIndex = 1, `numOrders'
* Construct `cppOrder' for use in the function and file names.
* {
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
* }

* define c++ preprocessor macros for the number of integration variables in each integrand
  #write <sector_`sectorID'.hpp> "#define sector_`sectorID'_order_`cppOrder'_numIV `numOccurringIVOrder`shiftedOrderIndex''#@SecDecInternalNewline@#"

* include the headers for all orders in this sector
  #write <sector_`sectorID'.hpp> "#include \"sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"

  #write <sector_`sectorID'.hpp> "#@SecDecInternalNewline@#"

* include contour deformation header (if needed)
  #If `contourDeformation'
    #write <sector_`sectorID'.hpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #EndIf

#EndDo

* open c++ namespace
#write <sector_`sectorID'.hpp> "namespace `name'#@SecDecInternalNewline@#"
#write <sector_`sectorID'.hpp> "{#@SecDecInternalNewline@#"

* define the data type of the container
#write <sector_`sectorID'.hpp> "`integrandContainerType'#@SecDecInternalNewline@#"

* write the c++ name of the conatainer instance
#write <sector_`sectorID'.hpp> "integrand_of_sector_`sectorID'#@SecDecInternalNewline@#"

* write constructor of the integrand container class
#write <sector_`sectorID'.hpp> "`integrandContainerInitializer';#@SecDecInternalNewline@#"

* close c++ namespace
#write <sector_`sectorID'.hpp> "};#@SecDecInternalNewline@#"

#write <sector_`sectorID'.hpp> "#@SecDecInternalNewline@#"

* undefine the c++ preprocessor macros for the number of integration variables
*{
#Do shiftedOrderIndex = 1, `numOrders'
* Construct `cppOrder' for use in the function and file names.
* {
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
* }

* undefine the c++ preprocessor macros
  #write <sector_`sectorID'.hpp> "#undef sector_`sectorID'_order_`cppOrder'_numIV#@SecDecInternalNewline@#"
#EndDo
*}

* finalize include guard
#write <sector_`sectorID'.hpp> "#endif#@SecDecInternalNewline@#"

.end
