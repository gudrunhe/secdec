#-
Off statistics;

* Define two general procedures that write c++ code to define and undefine
* c++ preprocessor varibables accessing a c++ array.
#procedure cppDefine(?FORMNames,cppArrayName,filename)
  #$counter = 0;
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <`filename'> "#define `varname' `cppArrayName'[`$counter']#@SecDecInternalNewline@#"
      #$counter = $counter + 1;
    #EndIf
  #EndDo
#endProcedure
#procedure cppUndefine(?FORMNames,filename)
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <`filename'> "#undef `varname'#@SecDecInternalNewline@#"
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
  #write <sector_`sectorID'_`cppOrder'.hpp> "#include \"functions.hpp\"#@SecDecInternalNewline@#"
  #If `contourDeformation'
    #write <sector_`sectorID'_`cppOrder'.hpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_complex_math.h>#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_linalg.h>#@SecDecInternalNewline@#"
  #EndIf
  #write <sector_`sectorID'_`cppOrder'.hpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "{#@SecDecInternalNewline@#"
  #If `contourDeformation'
    #write <sector_`sectorID'_`cppOrder'.hpp> "    secdecutil::SectorContainerWithDeformation<real_t, complex_t>::DeformedIntegrandFunction#@SecDecInternalNewline@#"
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
  #write <sector_`sectorID'_`cppOrder'.cpp> "    real_t const * const integration_variables,#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    real_t const * const real_parameters,#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "    complex_t const * const complex_parameters#@SecDecInternalNewline@#"
  #If `contourDeformation'
    #write <sector_`sectorID'_`cppOrder'.cpp> "    ,real_t const * const deformation_parameters,#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "    const real_t deformation_offset#@SecDecInternalNewline@#"
  #EndIf
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

* insert calI
  #call insertCalI

* insert the derivatives of the contour deformation Jacobian (should all be zero)
  #If `contourDeformation'
    #call insertContourdefJacobianDerivatives
  #EndIf
  .sort

* Replace calls to the deformation's Jacobian determinant "SecDecInternalContourdefJacobian" by symbols.
* {

  #If `contourDeformation'

    Local toOptimize = SecDecInternalsDUMMYtoOptimize;

    #redefine function "SecDecInternalContourdefJacobian"
    #$labelCounter = 0;

*   No need for "#Do depth = 0, `insertionDepth'" because the Jacobian determinant should only occur at top level.

*   Since we need intermediate ".sort" instructions, we cannot use the
*   "repeat" environment.
*   The following construction is suggested in the FORM documentation.

    #Do i = 1,1
*     set dollar variable
      .sort
      hide toOptimize;
      .sort
      if ( match(`function'(?SecDecInternalsDUMMY$args)) ) redefine i "0";
      .sort
      unhide toOptimize;
      .sort

*     The following "#if" evaluates to true only if there is still something to do.
      #If `i' == 0

        #$labelCounter = $labelCounter + 1;

        #Do replaceDepth = 0, `insertionDepth'
          #call beginArgumentDepth(`replaceDepth')
            Id `function'(`$args') = SecDecInternal`function'Call`$labelCounter';
          #call endArgumentDepth(`replaceDepth')
        #EndDo

        Id SecDecInternalsDUMMYtoOptimize = SecDecInternalsDUMMYtoOptimize + `function'($args) * SecDecInternalLabelSecDecInternalContourdefJacobianCall^`$labelCounter';

*       Remember how many integration variables are nonzero.
        #$index = 0;
        #$nonzeroCounter = 0;
        #redefine SecDecInternalContourdefLabel`$labelCounter'NonzeroIndices ""
        #Do arg = {`$args',}
          #If x`arg' != x
            #$index = $index + 1;
            #If `arg' != 0
              #$nonzeroCounter = $nonzeroCounter + 1;
              #If `$nonzeroCounter' != 1
                #redefine SecDecInternalContourdefLabel`$labelCounter'NonzeroIndices "`SecDecInternalContourdefLabel`$labelCounter'NonzeroIndices',"
              #EndIf
              #redefine SecDecInternalContourdefLabel`$labelCounter'NonzeroIndices "`SecDecInternalContourdefLabel`$labelCounter'NonzeroIndices'`$index'"
            #EndIf
          #EndIf
        #EndDo
        #redefine SecDecInternalContourdefLabel`$labelCounter'Dimensionality "`$nonzeroCounter'"

      #EndIf
      .sort
    #EndDo

    #redefine largestLabel`function' "`$labelCounter'"

*   insert the Jacobian matrix
    #call insertContourdefJacobianMatrix

  #EndIf

* }

* Replace calls to the deformation of the integration variables "SecDecInternalDeformed..." by symbols.
* {

  #If `contourDeformation'

    #Do depth = 0, `insertionDepth'
      #call beginArgumentDepth(`depth')

*       Do not expand functions to higher powers. --> Wrap into function "SecDecInternalIntPow"
*       example: "U(x,y,z)^2" --> "SecDecInternalIntPow(U(x,y,z),2)"
        repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) * SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) =
            SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), 2);
        repeat Id SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1?)
                * SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo2?)
                = SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1 + SecDecInternalsDUMMYexpo2);
        repeat Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 0) = 1;
        repeat Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
        repeat Id SecDecInternalIntPow(0, SecDecInternalsDUMMYexponent?) = 0;

*       Cancel ratios of functions and wrap denominators into the function "SecDecInternalDenominator".
*       example: "U(x,y,z)/U(x,y,z)^2" --> "SecDecInternalDenominator(U(x,y,z))"
        Denominators SecDecInternalDenominator;
        factarg,(-1),SecDecInternalDenominator;
        chainout SecDecInternalDenominator;
        repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;

      #call endArgumentDepth(`depth')
      .sort

    #EndDo

    #$IVindex = 0;
    #Do IV = {`integrationVariables'}
      #redefine function "SecDecInternalDeformed`IV'"
      #$labelCounter = 0;
      #$IVindex = $IVindex + 1;

*     Define a preprocessor variable to remove the deformed variable
*     if the undeformed equivalent is zero.
      #redefine matchArg ""
      #$counter = 0;
      #Do otherIV = {`integrationVariables'}
        #$counter = $counter + 1;
        #If `$counter' != 1
          #redefine matchArg "`matchArg',"
        #EndIf
        #If `$counter' == `$IVindex'
          #redefine matchArg "`matchArg'0"
        #Else
          #redefine matchArg "`matchArg'`otherIV'?"
        #EndIf
      #EndDo

*     remove deformed version if undeformed variable is set to zero
      #Do depth = 0, `insertionDepth'
        #call beginArgumentDepth(`depth')
          Id `function'(`matchArg') = 0;
        #call endArgumentDepth(`depth')
      #EndDo

      #Do depth = 0, `insertionDepth'

*       Since we need intermediate ".sort" instructions, we cannot use the
*       "repeat" environment.
*       The following construction is suggested in the FORM documentation.

        #Do i = 1,1
*         set dollar variable
          #call beginArgumentDepth(`depth')
            if ( match(`function'(?SecDecInternalsDUMMY$args)) ) redefine i "0";
          #call endArgumentDepth(`depth')
          .sort

*         The following "#if" evaluates to true only if there is still something to do.
          #If `i' == 0

            #$labelCounter = $labelCounter + 1;

            #Do replaceDepth = 0, `insertionDepth'
              #call beginArgumentDepth(`replaceDepth')
                Id `function'(`$args') = SecDecInternal`function'Call`$labelCounter';
              #call endArgumentDepth(`replaceDepth')
            #EndDo

            Id SecDecInternalsDUMMYtoOptimize = SecDecInternalsDUMMYtoOptimize +
                SecDecInternalLabel`function' ^ $labelCounter * `function'($args);

            .sort
          #EndIf
        #EndDo
      #EndDo

      #redefine largestLabel`function' "`$labelCounter'"

    #EndDo

    #call insertDeformedIntegrationVariables;
    .sort

  #EndIf

* }

* Explicitly insert the other functions defined in python.
* {

  #Do depth = 0, `insertionDepth'
    #call beginArgumentDepth(`depth')

*     Do not expand functions to higher powers. --> Wrap into function "SecDecInternalIntPow"
*     example: "U(x,y,z)^2" --> "SecDecInternalIntPow(U(x,y,z),2)"
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) * SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) =
          SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), 2);
        repeat Id SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1?)
                * SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo2?)
                = SecDecInternalIntPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1 + SecDecInternalsDUMMYexpo2);
      repeat Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 0) = 1;
      repeat Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
      repeat Id SecDecInternalIntPow(0, SecDecInternalsDUMMYexponent?) = 0;

*     Cancel ratios of functions and wrap denominators into the function "SecDecInternalDenominator".
*     example: "U(x,y,z)/U(x,y,z)^2" --> "SecDecInternalDenominator(U(x,y,z))"
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;

    #call endArgumentDepth(`depth')

    .sort

    #call beginArgumentDepth(`depth')
      #call insertOther
      #If `contourDeformation'
        #call insertDeformedIntegrationVariables
      #EndIf
    #call endArgumentDepth(`depth')
    .sort

*   some simplifications
    #call beginArgumentDepth(`depth')
      repeat Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 0) = 1;
      repeat Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
      repeat Id SecDecInternalIntPow(0, SecDecInternalsDUMMYexponent?) = 0;
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      repeat Id log(1) = 0;
      #If `contourDeformation'
        #Do function = {SecDecInternalExpMinusMuOverX,dSecDecInternalXExpMinusMuOverXd1,SecDecInternalXExpMinusMuOverX}
          repeat Id `function'(SecDecInternalsDUMMY?, 0) = 0;
        #EndDo
      #EndIf
      repeat Id SecDecInternalsDUMMY1? ^ SecDecInternalsDUMMY2?neg_ = SecDecInternalDenominator(SecDecInternalsDUMMY1) ^ (-SecDecInternalsDUMMY2);
      repeat Id 1/SecDecInternalsDUMMY? = SecDecInternalDenominator(SecDecInternalsDUMMY);
      repeat Id SecDecInternalsDUMMY? * SecDecInternalDenominator(SecDecInternalsDUMMY?) = 1;
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;
      repeat Id SecDecInternalDenominator(SecDecInternalsDUMMY?number_) = 1/SecDecInternalsDUMMY;
      #If `contourDeformation'
        repeat Id SecDecInternalRealPart(SecDecInternalsDUMMY?number_) = SecDecInternalsDUMMY;
      #EndIf
    #call endArgumentDepth(`depth')
    .sort
  #EndDo

* }

* translate sympy's imaginary unit to FORM's imaginary unit
  multiply replace_(I,i_);
  .sort

* simplify again
  #Do depth = 0, `insertionDepth'
    #call beginArgumentDepth(`depth')
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      repeat Id log(1) = 0;
      #If `contourDeformation'
        #Do function = {SecDecInternalExpMinusMuOverX,dSecDecInternalXExpMinusMuOverXd1,SecDecInternalXExpMinusMuOverX}
          repeat Id `function'(SecDecInternalsDUMMY?, 0) = 0;
        #EndDo
      #EndIf
      repeat Id SecDecInternalsDUMMY1? ^ SecDecInternalsDUMMY2?neg_ = SecDecInternalDenominator(SecDecInternalsDUMMY1) ^ (-SecDecInternalsDUMMY2);
      repeat Id 1/SecDecInternalsDUMMY? = SecDecInternalDenominator(SecDecInternalsDUMMY);
      repeat Id SecDecInternalsDUMMY? * SecDecInternalDenominator(SecDecInternalsDUMMY?) = 1;
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;
      repeat Id SecDecInternalDenominator(SecDecInternalsDUMMY?number_) = 1/SecDecInternalsDUMMY;
      #If `contourDeformation'
        repeat Id SecDecInternalRealPart(SecDecInternalsDUMMY?number_) = SecDecInternalsDUMMY;
      #EndIf
    #call endArgumentDepth(`depth')
  #EndDo
  repeat;
    #Do depth = 0, `insertionDepth'
      #call beginArgumentDepth(`depth')
        Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 0) = 1;
        Id SecDecInternalIntPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
        Id SecDecInternalIntPow(0, SecDecInternalsDUMMYexponent?) = 0;
      #call endArgumentDepth(`depth')
    #EndDo
  endRepeat;
  .sort

* Replace all function calls by symbols for simultaneous optimization.
* {

  #If `contourDeformation' == 0
    Local toOptimize = SecDecInternalsDUMMYtoOptimize;
  #EndIf

  #redefine functionsToReplace "`functions',log,SecDecInternalIntPow,SecDecInternalDenominator"
  #If `contourDeformation'
    #redefine functionsToReplace "SecDecInternalRealPart,`functionsToReplace'"
  #EndIf

  #Do function = {`functionsToReplace'}
    #redefine largestLabel`function' "0"
  #EndDo

  #Do j = 1,1

    #Do function = {`functionsToReplace'}
      #$labelCounter = `largestLabel`function'';

      #Do depth = 0, `insertionDepth'

*       Since we need intermediate ".sort" instructions, we cannot use the
*       "repeat" environment.
*       The following construction is suggested in the FORM documentation.

        #Do i = 1,1
*         set dollar variable
          #call beginArgumentDepth(`depth')
            if ( match(`function'(?SecDecInternalsDUMMY$args)) ) redefine i "0";
          #call endArgumentDepth(`depth')
          .sort

*         The following "#if" evaluates to true only if there is still something to do.
          #If `i' == 0

            #redefine j "0"

            #$labelCounter = $labelCounter + 1;

            L arguments = SecDecInternalfDUMMYarguments(`$args');

            #Do replaceDepth = 0, `insertionDepth'
              #call beginArgumentDepth(`replaceDepth')
                Id `function'(`$args') = SecDecInternal`function'Call`$labelCounter';
              #call endArgumentDepth(`replaceDepth')
            #EndDo

            repeat Id SecDecInternalfDUMMYarguments(SecDecInternalsDUMMY?, ?otherArgs) =
                SecDecInternalLabel`function'Call`$labelCounter'Arg * (SecDecInternalsDUMMY + SecDecInternalfDUMMYarguments(?otherArgs));

*           Define `$argCounter' by loking at the term with the empty function "SecDecInternalfDUMMYarguments"
            Id SecDecInternalfDUMMYarguments * SecDecInternalLabel`function'Call`$labelCounter'Arg ^ SecDecInternalsDUMMYexponent?$argCounter = 0;
            .sort

*           Add all arguments to top level polynomial for simultaneous optimization.
            Id SecDecInternalsDUMMYtoOptimize = SecDecInternalsDUMMYtoOptimize + arguments;

            #redefine numberOfArgs`function'Label`$labelCounter' "`$argCounter'"

            .sort
          #EndIf
        #EndDo
      #EndDo

      #redefine largestLabel`function' "`$labelCounter'"

    #EndDo

  #EndDo

  Id SecDecInternalsDUMMYtoOptimize = expression;
  .sort

  drop expression, arguments;
  .sort
* }

* Find and count the occurring integration variables.
* {
  hide; nhide toOptimize;
  .sort

  #redefine occurringIntegrationVariables ""
  #redefine occurringIntegrationVariableIndices ""
  #redefine absentIntegrationVariables ""
  #redefine absentIntegrationVariableIndices ""
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

  unhide;
  .sort
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
  AntiBracket `integrationVariables', `realParameters', `complexParameters'
  #If `contourDeformation'
    , `deformationParameters', SecDecInternalMu
  #EndIf
  ;
  Format O`optimizationLevel';
  .sort
  #optimize toOptimize

* Define the integration variables and parameters as c preprocessor variables
* (The integrand function in c takes them packed into an array).
* {
* "Format rational": Need the indices as integers.
  Format rational;

* call the general procedure to write the corresponding c++ code define in the beginning of this file
  #call cppDefine(`occurringIntegrationVariables',integration_variables,sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`realParameters',real_parameters,sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`complexParameters',complex_parameters,sector_`sectorID'_`cppOrder'.cpp)
  #If `contourDeformation'
    #call cppDefine(`occurringDeformationParameters',deformation_parameters,sector_`sectorID'_`cppOrder'.cpp)
    #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalMu deformation_offset#@SecDecInternalNewline@#"
  #EndIf
* }

* Processing denominators in FORM is easiest if packed into a function.
* Define that function as c preprocessor macro.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalDenominator(x) 1./(x)#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) std::real(x)#@SecDecInternalNewline@#"

* Define "SecDecInternalAbbreviation[0]" as c preprocessor variable "tmp".
* Since FORM does not use "SecDecInternalAbbreviation[0]", we can use it.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* write Abbreviations in c format
  Format float 20;
  Format C;
  #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* c++ define the calls to "SecDecInternalRealPart" BEFORE computing the Jacobian determinant
* {
  #If `contourDeformation'

    #redefine function "SecDecInternalRealPart"

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
        #write <sector_`sectorID'_`cppOrder'.cpp> "%%E"  arg`argIndex'(#@no_split_expression@#)
        #If `argIndex' != `numberOfArgs`function'Label`callIndex''
          #write <sector_`sectorID'_`cppOrder'.cpp> ","
        #EndIf
      #EndDo
      #write <sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
      multiply replace_(SecDecInternalLabel`function'Call`callIndex'Arg,0);
      .sort
      #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
        drop arg`argIndex';
      #EndDo
    #EndDo
    .sort

  #EndIf
* }

* Define the calls to the contour deformation "SecDecInternalDeformed..."
* and "SecDecInternalContourdefJacobian".
* {

  #If `contourDeformation'

*   c++ define the defomed integration variables
    #Do IV = {`occurringIntegrationVariables',}
      #If x`IV' != x
        #Do callIndex = 1, `largestLabelSecDecInternalDeformed`IV''

          Bracket SecDecInternalLabelSecDecInternalDeformed`IV';
          .sort
          L deformedIV = toOptimize[SecDecInternalLabelSecDecInternalDeformed`IV' ^ `callIndex'];
          .sort
          #write <sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalSecDecInternalDeformed`IV'Call`callIndex' = %%e" deformedIV(#@no_split_expression@#)
          #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

        #EndDo

*       Remove label once all deformations of this variable are parsed.
        multiply replace_(SecDecInternalLabelSecDecInternalDeformed`IV', 0);
        .sort

      #EndIf

    drop deformedIV;
    #EndDo

    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

    Format float 20;
    Format C;

*   c++ define the Jacobian determinant of the contour deformation.
*   We use th gsl to compute it numerically out of the Jacobian matrix.
*   {
    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "// variables required to compute a complex numerical determinant using the gsl#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex *Jacobian_matrix;#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_complex Jacobian_determinant;#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_permutation * Jacobian_permutation;#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "int Jacobian_signum;#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "unsigned Jacobian_dimensionality;#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

    #Do callIndex = 1, `largestLabelSecDecInternalContourdefJacobian'
      #write <sector_`sectorID'_`cppOrder'.cpp> "// begin code for numerical Jacobian determinant `callIndex'#@SecDecInternalNewline@#"
      #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

      #If `SecDecInternalContourdefLabel`callIndex'Dimensionality' > 0

*       allocate memory
        #write <sector_`sectorID'_`cppOrder'.cpp> "Jacobian_dimensionality = `SecDecInternalContourdefLabel`callIndex'Dimensionality';#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "Jacobian_matrix = gsl_matrix_complex_alloc("
        #write <sector_`sectorID'_`cppOrder'.cpp> "Jacobian_dimensionality,Jacobian_dimensionality"
        #write <sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "Jacobian_permutation = gsl_permutation_alloc(Jacobian_dimensionality);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

*       fill the gsl matrix
        #$i = -1;
        #Do idx1 = {`SecDecInternalContourdefLabel`callIndex'NonzeroIndices',}
          #If x`idx1' != x
            #$i = $i + 1;
            #$j = -1;
            #Do idx2 = {`SecDecInternalContourdefLabel`callIndex'NonzeroIndices',}
              #If x`idx2' != x
                Format float 20;
                Format C;
                #$j = $j + 1;
                Bracket SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ, SecDecInternalLabelSecDecInternalContourdefJacobianCall;
                .sort
                L expr = toOptimize[SecDecInternalLabelJacobianMatrixI^`idx1' * SecDecInternalLabelJacobianMatrixJ^`idx2' * SecDecInternalLabelSecDecInternalContourdefJacobianCall^`callIndex'];
                .sort
                #write <sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e" expr(#@no_split_expression@#)
                #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex_set#@SecDecInternalNewline@#"
                #write <sector_`sectorID'_`cppOrder'.cpp> "(#@SecDecInternalNewline@#"
                #write <sector_`sectorID'_`cppOrder'.cpp> "    Jacobian_matrix,#@SecDecInternalNewline@#"
                Format rational;
                #write <sector_`sectorID'_`cppOrder'.cpp> "    `$i', `$j',#@SecDecInternalNewline@#"
                Format float 20;
                Format C;
                #write <sector_`sectorID'_`cppOrder'.cpp> "    {tmp.real(), tmp.imag()}#@SecDecInternalNewline@#"
                #write <sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
              #EndIf
            #EndDo
          #EndIf
        #EndDo

*       calculate the determinant numerically using the gsl
        #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_linalg_complex_LU_decomp(Jacobian_matrix, Jacobian_permutation, &Jacobian_signum);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "Jacobian_determinant = gsl_linalg_complex_LU_det(Jacobian_matrix, Jacobian_signum);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalSecDecInternalContourdefJacobianCall`callIndex' ="
        #write <sector_`sectorID'_`cppOrder'.cpp> "{GSL_REAL(Jacobian_determinant), GSL_IMAG(Jacobian_determinant)};#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

*       free manually allocated memory
        #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_permutation_free(Jacobian_permutation);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex_free(Jacobian_matrix);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

      #Else

*       dimensionality is zero --> no transformation --> determinant is one
        #write <sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalSecDecInternalContourdefJacobianCall`callIndex' = 1;#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

      #EndIf

      #write <sector_`sectorID'_`cppOrder'.cpp> "// end code for numerical Jacobian determinant `callIndex'#@SecDecInternalNewline@#"
    #EndDo

    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
*   }

*   remove the parsed labels
    Id SecDecInternalLabelSecDecInternalContourdefJacobianCall = 0;
    .sort
    drop expr;

  #EndIf

* }

* Define the function calls replaced by symbols.
* The difficulty here is that the function arguments may
* refer to other function calls. --> We must order the
* definitions in the c++ file properly.
* {

* Keep track of function calls that are not written to
* the c++ file yet.
  L unparsed = SecDecInternalsDUMMYUnparsedAppendix;
  #Do function = {`functions',log,SecDecInternalIntPow,SecDecInternalDenominator}
    #Do callIndex = 1, `largestLabel`function''
      Id SecDecInternalsDUMMYUnparsedAppendix = SecDecInternalsDUMMYUnparsedAppendix + SecDecInternal`function'Call`callIndex'Unparsed;
    #EndDo
  #EndDo

  #Do i = 1,1
    #Do function = {`functions',log,SecDecInternalIntPow,SecDecInternalDenominator}
      #Do callIndex = 1, `largestLabel`function''
        B SecDecInternalLabel`function'Call`callIndex'Arg;
        .sort
        #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
          L arg`argIndex' = toOptimize[SecDecInternalLabel`function'Call`callIndex'Arg ^ `argIndex'];
        #EndDo
        .sort

*       We do not want to define any call more than once.
        #redefine alreadyParsed "1"
        if ( occurs(SecDecInternal`function'Call`callIndex'Unparsed) ) redefine alreadyParsed "0";
        .sort
        #If `alreadyParsed' == 0

*         We can write the call under consideration to the c++ file only
*         if all dependent calls are already written.
          #redefine dependenciesDone "1"
          hide; nhide unparsed;
          #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
            nhide arg`argIndex';
          #EndDo
          .sort
          #Do innerFunction = {`functions',log,SecDecInternalIntPow,SecDecInternalDenominator}
            #Do innerCallIndex = 1, `largestLabel`innerFunction''
              #redefine dependsOnInnerFunctionCall "0"
              if ( occurs(SecDecInternal`innerFunction'Call`innerCallIndex') ) redefine dependsOnInnerFunctionCall "1";
              .sort
              #If `dependsOnInnerFunctionCall'
                if ( occurs(SecDecInternal`innerFunction'Call`innerCallIndex'Unparsed) );
                  redefine dependenciesDone "0";
                  redefine i "0";
                endif;
              #EndIf
            #EndDo
          #EndDo
          .sort
          unhide;
          .sort

          #If `dependenciesDone'
            #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternal`function'Call`callIndex' = "
            #write <sector_`sectorID'_`cppOrder'.cpp> "`function'("
            #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
              #write <sector_`sectorID'_`cppOrder'.cpp> "%%E"  arg`argIndex'(#@no_split_expression@#)
              #If `argIndex' != `numberOfArgs`function'Label`callIndex''
                #write <sector_`sectorID'_`cppOrder'.cpp> ","
              #EndIf
            #EndDo
            #write <sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
            multiply replace_(SecDecInternal`function'Call`callIndex'Unparsed,0 , SecDecInternalLabel`function'Call`callIndex'Arg,0 );
            .sort
          #EndIf
          #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
            drop arg`argIndex';
          #EndDo
        #EndIf
      #EndDo
    #EndDo
  #EndDo

  .sort
  drop unparsed;
* }

* write the integrand
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "return %%e#@SecDecInternalNewline@#" toOptimize(#@no_split_expression@#)
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables',sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`realParameters',sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`complexParameters',sector_`sectorID'_`cppOrder'.cpp)
  #If `contourDeformation'
    #call cppUndefine(`occurringDeformationParameters',sector_`sectorID'_`cppOrder'.cpp)
    #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalMu#@SecDecInternalNewline@#"
  #EndIf
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalDenominator#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the c++ function and namespaces
  #write <sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* write the contour deformation optimize functions if required
  #If `contourDeformation'
    #include write_contour_deformation.frm
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

* include contour deformation and optimize deformation_parameter headers (if needed)
  #If `contourDeformation'
    #write <sector_`sectorID'.hpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
    #write <sector_`sectorID'.hpp> "#include \"optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #EndIf

  #write <sector_`sectorID'.hpp> "#@SecDecInternalNewline@#"

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
