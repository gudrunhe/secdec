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

* Define a simplification procedure
#procedure simplify
  #Do depth = 0, `insertionDepth'
    #call beginArgumentDepth(`depth')
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      repeat Id log(1) = 0;
      repeat Id SecDecInternalsDUMMY? * SecDecInternalsDUMMY? = SecDecInternalPow(SecDecInternalsDUMMY, 2);
      repeat Id SecDecInternalsDUMMYbase? * SecDecInternalPow(SecDecInternalsDUMMYbase?, SecDecInternalsDUMMYexponent?) =
        SecDecInternalPow(SecDecInternalsDUMMYbase, SecDecInternalsDUMMYexponent + 1);
      repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, SecDecInternalsDUMMYexponent1?) * SecDecInternalPow(SecDecInternalsDUMMYbase?, SecDecInternalsDUMMYexponent2?) =
        SecDecInternalPow(SecDecInternalsDUMMYbase, SecDecInternalsDUMMYexponent1 + SecDecInternalsDUMMYexponent2);
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
        Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 0) = 1;
        Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
        Id SecDecInternalPow(0, SecDecInternalsDUMMYexponent?) = 0;
        Id SecDecInternalPow(1, SecDecInternalsDUMMYexponent?) = 1;
      #call endArgumentDepth(`depth')
    #EndDo
  endRepeat;
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
    #write <sector_`sectorID'_`cppOrder'.cpp> "    ,real_t const * const deformation_parameters#@SecDecInternalNewline@#"
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

* Find the calls to the contour deformation polynomials that need a sign check.
* The call differ only in which Feynman paramters are set to zero. We investigate
* that by looking at the calls to "SecDecInternalCalI".
* {
  #If `contourDeformation'

    Local tmp = expression;
    Local signCheck = SecDecInternalsDUMMYsignCheck;
    .sort
    hide; nhide tmp, signCheck;
    .sort

    #$labelCounter = 0;

    #Do i = 1,1
*     set dollar variable
      .sort
      hide signCheck;
      .sort
      if ( match(SecDecInternalCalI(?SecDecInternalsDUMMY$args)) ) redefine i "0";
      .sort
      unhide signCheck;
      .sort

*     The following "#if" evaluates to true only if there is still something to do.
      #If `i' == 0

        #$labelCounter = $labelCounter + 1;

        Id SecDecInternalCalI(`$args') = 0;

*       The sign check to be performed is:
*       "SecDecInternalContourDeformationPolynomial(<deformed integration variables>) - SecDecInternalContourDeformationPolynomial(<undeformed integration variables>) <= 0"
*       {
        Id SecDecInternalsDUMMYsignCheck = SecDecInternalsDUMMYsignCheck +
            SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheck`$labelCounter' *
            (
              SecDecInternalfDUMMYdeformedContourDeformationPolynomial($args) - `SecDecInternalContourDeformationPolynomial'($args)
            );

*       Define a variables that only keeps the integration variables in "$args"
*       {
        #redefine argsWithoutRegulators ""
        #$counter = 0;
        #Do arg = {`$args'}
          #$counter = $counter + 1;
          #If `$counter' <= `numIV'
            #If `$counter' != 1
              #redefine argsWithoutRegulators "`argsWithoutRegulators',"
            #EndIf
            #redefine argsWithoutRegulators "`argsWithoutRegulators'`arg'"
          #EndIf
        #EndDo
*       }

        #Do IV = {`integrationVariables'}
          Id SecDecInternalfDUMMYdeformedContourDeformationPolynomial(?frontArgs,`IV',?backArgs) =
             SecDecInternalfDUMMYdeformedContourDeformationPolynomial(?frontArgs,SecDecInternalDeformed`IV'(`argsWithoutRegulators'),?backArgs);
        #EndDo

        Id SecDecInternalfDUMMYdeformedContourDeformationPolynomial(?args) = `SecDecInternalContourDeformationPolynomial'(?args);

        .sort

*       }

      #EndIf
    #EndDo

    #redefine numberOfRequiredSignChecks "`$labelCounter'"

    unhide;
    drop tmp;
    multiply replace_(SecDecInternalsDUMMYsignCheck,0);
    .sort

  #EndIf
* }

* insert calI
  #call insertCalI
  .sort

  Local toOptimize = SecDecInternalsDUMMYtoOptimize;

* Replace calls to the deformation's Jacobian determinant "SecDecInternalContourdefJacobian" by symbols.
* {

  #If `contourDeformation'

    #Do function = {`contourdefJacobianFunctions'}
      #$labelCounter = 0;

*     No need for "#Do depth = 0, `insertionDepth'" because the Jacobian determinant should only occur at top level.

*     Since we need intermediate ".sort" instructions, we cannot use the
*     "repeat" environment.
*     The following construction is suggested in the FORM documentation.

      #Do i = 1,1
*       set dollar variable
        .sort
        hide toOptimize;
        .sort
        if ( match(`function'(?SecDecInternalsDUMMY$args)) ) redefine i "0";
        .sort
        unhide toOptimize;
        .sort

*       The following "#if" evaluates to true only if there is still something to do.
        #If `i' == 0

          #$labelCounter = $labelCounter + 1;

          Id `function'(`$args') = SecDecInternal`function'Call`$labelCounter';

          Id SecDecInternalsDUMMYtoOptimize = SecDecInternalsDUMMYtoOptimize +
              SecDecInternalLabel`function'Call`$labelCounter'Arg * `function'($args);

*           no arguments after "insertContourdefJacobianDerivatives"
            #redefine numberOfArgs`function'Label`$labelCounter' "0"

        #EndIf
        .sort
      #EndDo

      #redefine largestLabel`function' "`$labelCounter'"

*     insert the contour deformation Jacobian and its derivatives
      #call insertContourdefJacobianDerivatives
      .sort

    #EndDo

  #EndIf

* }

* Replace calls to the deformation of the integration variables "SecDecInternalDeformed..."
* (and appearing derivatives) by symbols.
* {

  #If `contourDeformation'

    #Do depth = 0, `insertionDepth'
      #call beginArgumentDepth(`depth')

*       Do not expand functions to higher powers. --> Wrap into function "SecDecInternalPow"
*       example: "U(x,y,z)^2" --> "SecDecInternalPow(U(x,y,z),2)"
        repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) * SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) =
            SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), 2);
        repeat Id SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1?)
                * SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo2?)
                = SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1 + SecDecInternalsDUMMYexpo2);
        repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 0) = 1;
        repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
        repeat Id SecDecInternalPow(0, SecDecInternalsDUMMYexponent?) = 0;

*       Wrap noninteger powers into the function pow.
        repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) ^ SecDecInternalsDUMMYExponent? =
            SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYExponent);

*       Cancel ratios of functions and wrap denominators into the function "SecDecInternalDenominator".
*       example: "U(x,y,z)/U(x,y,z)^2" --> "SecDecInternalDenominator(U(x,y,z))"
        Denominators SecDecInternalDenominator;
        factarg,(-1),SecDecInternalDenominator;
        chainout SecDecInternalDenominator;
        repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;

      #call endArgumentDepth(`depth')
      .sort

    #EndDo

    #Do function = {`deformedIntegrationVariableDerivativeFunctions'}
      #$labelCounter = 0;
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

*           check if the requested call is zero
*           {
            Format rational;
            .sort
            Local zeroCheck = `function'($args) * SecDecInternalsDUMMYZeroCheck;
            .sort
            hide; nhide zeroCheck;
            .sort

            #call insertDeformedIntegrationVariables
            Argument SecDecInternalRealPart;
              #call insertOther
            EndArgument;
            Id SecDecInternalRealPart(SecDecInternalsDUMMY?number_) = SecDecInternalsDUMMY;
*           Remove the deformation parameters of absent integration variables.
*           {
            #$argCounter = 0;
            #Do arg = {`$args',}
              #If x`arg' != x
                #If `arg' == 0
                  multiply replace_(SecDecInternalLambda`$argCounter',0);
                  .sort
                #EndIf
                #$argCounter = $argCounter + 1;
              #EndIf
            #EndDo
*           }
            multiply replace_(I,i_);
            .sort

            #redefine callIsZero "1"
            if ( occurs(SecDecInternalsDUMMYZeroCheck) ) redefine callIsZero "0";
            .sort

            unhide;
            drop zeroCheck;
            .sort
*           }

            #If `callIsZero'

              multiply replace_(SecDecInternal`function'Call`$labelCounter',0);
              #$labelCounter = $labelCounter - 1;

            #Else

              Id SecDecInternalsDUMMYtoOptimize = SecDecInternalsDUMMYtoOptimize +
                  SecDecInternalLabel`function'Call`$labelCounter'Arg * SecDecInternalfDUMMY`function'($args);

*             no arguments after "insertDeformedIntegrationVariables"
              #redefine numberOfArgs`function'Label`$labelCounter' "0"

            #EndIf

            .sort

          #EndIf
        #EndDo

      repeat Id SecDecInternalfDUMMY`function'(?arguments) = `function'(?arguments);

      .sort

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

*     Do not expand functions to higher powers. --> Wrap into function "SecDecInternalPow"
*     example: "U(x,y,z)^2" --> "SecDecInternalPow(U(x,y,z),2)"
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) * SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) =
          SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), 2);
        repeat Id SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1?)
                * SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo2?)
                = SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYexpo1 + SecDecInternalsDUMMYexpo2);
      repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 0) = 1;
      repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
      repeat Id SecDecInternalPow(0, SecDecInternalsDUMMYexponent?) = 0;

*     Wrap noninteger powers into the function pow.
      repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMYArgs) ^ SecDecInternalsDUMMYExponent? =
          SecDecInternalPow(SecDecInternalfDUMMY(?SecDecInternalsDUMMYArgs), SecDecInternalsDUMMYExponent);

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
      repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 0) = 1;
      repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
      repeat Id SecDecInternalPow(0, SecDecInternalsDUMMYexponent?) = 0;
      Denominators SecDecInternalDenominator;
      factarg,(-1),SecDecInternalDenominator;
      chainout SecDecInternalDenominator;
      repeat Id log(1) = 0;
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

  #call simplify
  .sort

* Find and count the occurring integration variables.
* {
* "Format rational" because we need the dollar variables as integers
  Format rational;
  hide; nhide toOptimize, expression;
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

    #redefine integrationVariable`$currentIVIndex' "`IV'"

  #EndDo

  #redefine numOccurringIVOrder`shiftedOrderIndex' "`$counterOccur'"
  #redefine integrationVariable "SecDecInternalError"

  unhide;
  .sort
* }

* Specify the occurring/absent deformation parameters.
* {
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

    #redefine absentDeformationParameters ""
    #$counter = 1;
    #Do var = {`absentIntegrationVariableIndices',}
      #If x`var' != x
        #If `$counter' != 1
           #redefine absentDeformationParameters "`absentDeformationParameters',"
        #Endif
        #redefine absentDeformationParameters "`absentDeformationParameters'SecDecInternalLambda`var'"
        #$counter = $counter + 1;
      #EndIf
    #EndDo
  #EndIf
* }

* Remove the deformation parameters of absent integration variables.
* {
  #If `contourDeformation'
    #redefine replaceArg ""
    #$counter = 1;
    #Do var = {`absentDeformationParameters',}
      #If x`var' != x
        #If `$counter' != 1
           #redefine replaceArg "`replaceArg' ,"
        #Endif
        #redefine replaceArg "`replaceArg' `var',0"
        #$counter = $counter + 1;
      #EndIf
    #EndDo
  multiply replace_(`replaceArg');
  .sort
  #EndIf
* }

  #call simplify
  .sort

* Check if any call to the deformation of the integration variables "SecDecInternalDeformed..."
* (and appearing derivatives) is zero or one and can therefore be removed.
* {

  #If `contourDeformation'

    #Do function = {`deformedIntegrationVariableDerivativeFunctions'}
      #Do callIndex = 1, `largestLabel`function''

        B SecDecInternalLabel`function'Call`callIndex'Arg;
        .sort

        L zeroCheck = toOptimize[SecDecInternalLabel`function'Call`callIndex'Arg] * SecDecInternalsDUMMYZeroCheck - SecDecInternalsDUMMYOneCheck;
        .sort
        hide; nhide zeroCheck;
        .sort

        #redefine callIsZero "1"
        #redefine callIsOne "1"
        if ( occurs(SecDecInternalsDUMMYZeroCheck) ) redefine callIsZero "0";
        .sort

        #If `callIsZero'
          #redefine callIsOne "0"
        #Else
          Id SecDecInternalsDUMMYZeroCheck = SecDecInternalsDUMMYOneCheck;
          .sort
          if ( occurs(SecDecInternalsDUMMYOneCheck) ) redefine callIsOne "0";
          .sort
        #EndIf

        unhide;
        drop zeroCheck;
        .sort

        #If `callIsZero'
          multiply replace_(SecDecInternal`function'Call`callIndex',0);
          .sort
        #ElseIf `callIsOne'
          #Do depth = 0, `insertionDepth'
            #call beginArgumentDepth(`depth')
              Id SecDecInternal`function'Call`callIndex' = 1;
            #call endArgumentDepth(`depth')
          #EndDo
          .sort
        #EndIf

      #EndDo
    #EndDo

  #EndIf

* }

  #call simplify
  .sort

* Replace all function calls by symbols for simultaneous optimization.
* {

  #redefine functionsToReplace "`functions',log,SecDecInternalPow,SecDecInternalDenominator"
  #If `contourDeformation'
    #redefine functionsToReplace "SecDecInternalRealPart,`functionsToReplace'"
  #EndIf

  #Do function = {`functionsToReplace'}
    #IfNDef `largestLabel`function''
      #redefine largestLabel`function' "0"
    #EndIf
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

  #If `contourDeformation'
    Id SecDecInternalsDUMMYtoOptimize = expression + signCheck;
  #Else
    Id SecDecInternalsDUMMYtoOptimize = expression;
  #EndIf
  .sort

  drop expression, arguments;
  #If `contourDeformation'
    drop signCheck;
  #EndIf
  .sort
* }

* Simultaneously optimize the integrand and all occurring function arguments.
  AntiBracket `integrationVariables', `realParameters', `complexParameters'
  #If `contourDeformation'
    , `deformationParameters'
  #EndIf
  ;
  Format O`optimizationLevel';
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviation;
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
  #EndIf
* }

* Processing denominators in FORM is easiest if packed into a function.
* Define that function as c preprocessor macro.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalDenominator(x) 1./(x)#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) std::real(x)#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalImagPart(x) std::imag(x)#@SecDecInternalNewline@#"

* Define "SecDecInternalAbbreviation[0]" as c preprocessor variable "tmp".
* Since FORM does not use "SecDecInternalAbbreviation[0]", we can use it.
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

* Define a function preprocessor macro to handle the abbrevition vectors
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalAbbreviations`shiftedOrderIndex'(ID) SecDecInternalAbbreviation[ID]#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* define the abbreviations in c
  Format float 20;
  Format C;
  #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternalSecondAbbreviation[sector_`sectorID'_order_`cppOrder'_optimmaxvar_second + 1];#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* write Abbreviations in c format
  #write <sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* Reload "toOptmize" such that the optimization symbols of the first
* optimization are promoted to regular FORM variables.
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviations`shiftedOrderIndex';
  #$toOptimize = toOptimize;
  #clearoptimize
  drop toOptimize;
  .sort
  AutoDeclare CFunctions SecDecInternalAbbreviations`shiftedOrderIndex';
  Format O0;
  Format normal;
  Local toOptimize = `$toOptimize';
  .sort

* Define the function calls replaced by symbols.
* The difficulty here is that the function arguments may
* refer to other function calls. --> We must order the
* definitions in the c++ file properly.
* {

* Keep track of function calls that are not written to
* the c++ file yet.

  L unparsed = SecDecInternalsDUMMYUnparsedAppendix;
  #redefine functionsToReplace "`functions',log,SecDecInternalPow,SecDecInternalDenominator"
  #If `contourDeformation'
    #redefine functionsToReplace "SecDecInternalRealPart,`contourdefJacobianFunctions',`deformedIntegrationVariableDerivativeFunctions',`functionsToReplace'"
  #EndIf

  #Do function = {`functionsToReplace'}
    #Do callIndex = 1, `largestLabel`function''
      Id SecDecInternalsDUMMYUnparsedAppendix = SecDecInternalsDUMMYUnparsedAppendix + SecDecInternal`function'Call`callIndex'Unparsed;
    #EndDo
  #EndDo

  #$globalLabelCounter = 1;
  #$globalIDCounter = 1;
  #redefine numberOfSecondAbbreviations "0"

  #Do i = 1,1
    Format normal;
    .sort
    #redefine parseNextGlobalIDMin "`$globalIDCounter'"
    #redefine parseNextGlobalIDMax "0"
    Local parseNext = SecDecInternalsDUMMYparseNext;
    #Do function = {`functionsToReplace'}
      #Do callIndex = 1, `largestLabel`function''
        B SecDecInternalLabel`function'Call`callIndex'Arg;
        .sort
        #If `numberOfArgs`function'Label`callIndex'' == 0
          L expr = toOptimize[SecDecInternalLabel`function'Call`callIndex'Arg];
        #Else
          #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
            L arg`argIndex' = toOptimize[SecDecInternalLabel`function'Call`callIndex'Arg ^ `argIndex'];
          #EndDo
        #EndIf
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
          #If `numberOfArgs`function'Label`callIndex'' == 0
            nhide expr;
          #Else
            #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
              nhide arg`argIndex';
            #EndDo
          #EndIf
          .sort
          #Do innerFunction = {`functionsToReplace'}
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

*           Define the variable only if it is used.
            #redefine needThisCall "0"
            if ( occurs(SecDecInternal`function'Call`callIndex') ) redefine needThisCall "1";
            .sort

*           Prepare for another optimization
*           {
            #If `needThisCall'

              #redefine parseNextGlobalIDMax "`$globalIDCounter'"
              #redefine parseNextFunction`$globalIDCounter' "`function'"
              #redefine parseNextCallIndex`$globalIDCounter' "`callIndex'"
              #$globalIDCounter = $globalIDCounter + 1;
              #redefine SecDecInternalLabel`function'Call`callIndex'GlobalIDMin "`$globalLabelCounter'"
              Format normal;

              #If `numberOfArgs`function'Label`callIndex'' == 0
                Id SecDecInternalsDUMMYparseNext = SecDecInternalsDUMMYparseNext +
                  expr * SecDecInternalLabelSecDecInternalGeneral ^ `$globalLabelCounter';
                #$globalLabelCounter = $globalLabelCounter + 1;
              #Else
                #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
                  Id SecDecInternalsDUMMYparseNext = SecDecInternalsDUMMYparseNext +
                    arg`argIndex' * SecDecInternalLabelSecDecInternalGeneral ^ `$globalLabelCounter';
                  #$globalLabelCounter = $globalLabelCounter + 1;
                #EndDo

              #EndIf

            #Else

              multiply replace_(SecDecInternal`function'Call`callIndex'Unparsed,0 , SecDecInternalLabel`function'Call`callIndex'Arg,0 );

            #EndIf
*           }

          #EndIf

          .sort
          #If `numberOfArgs`function'Label`callIndex'' == 0
            drop expr;
          #Else
            #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
              drop arg`argIndex';
            #EndDo
          #EndIf

        #EndIf
      #EndDo
    #EndDo

*   optimize and write to c file
*   {
    Format float 20;
    Format C;
    Format O`optimizationLevel';
    multiply replace_(SecDecInternalsDUMMYparseNext,0);
    Bracket SecDecInternalLabelSecDecInternalGeneral;
    .sort
    ExtraSymbols,array,SecDecInternalSecondAbbreviation;
    #optimize parseNext

    #If `optimmaxvar_' > `numberOfSecondAbbreviations'
      #redefine numberOfSecondAbbreviations "`optimmaxvar_'"
    #EndIf

    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@##@SecDecInternalNewline@#// begin next dependency level#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

    #Do globalID = `parseNextGlobalIDMin', `parseNextGlobalIDMax'
      #redefine function "`parseNextFunction`globalID''"
      #redefine callIndex "`parseNextCallIndex`globalID''"

      #If x`function' != x
        #write <sector_`sectorID'_`cppOrder'.cpp> "auto SecDecInternal`function'Call`callIndex' = "
        #If `numberOfArgs`function'Label`callIndex'' == 0
          Bracket SecDecInternalLabelSecDecInternalGeneral;
          .sort
          Local expr = parseNext[SecDecInternalLabelSecDecInternalGeneral ^ `SecDecInternalLabel`function'Call`callIndex'GlobalIDMin'];
          .sort
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%e#@SecDecInternalNewline@#"  expr(#@no_split_expression@#)
          drop expr;
        #ElseIf `function' == SecDecInternalPow
          Bracket SecDecInternalLabelSecDecInternalGeneral;
          .sort
          Local expr = parseNext[SecDecInternalLabelSecDecInternalGeneral ^ (`SecDecInternalLabel`function'Call`callIndex'GlobalIDMin')]
                     ^ parseNext[SecDecInternalLabelSecDecInternalGeneral ^ (`SecDecInternalLabel`function'Call`callIndex'GlobalIDMin'+1)];
          .sort
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%e#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
          drop expr;
        #Else
          #write <sector_`sectorID'_`cppOrder'.cpp> "`function'("
          #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
            Bracket SecDecInternalLabelSecDecInternalGeneral;
            .sort
            Local arg`argIndex' = parseNext[SecDecInternalLabelSecDecInternalGeneral ^ (`SecDecInternalLabel`function'Call`callIndex'GlobalIDMin'+`argIndex'-1)];
            .sort
            #write <sector_`sectorID'_`cppOrder'.cpp> "%%E"  arg`argIndex'(#@no_split_expression@#)
            drop arg`argIndex';
            #If `argIndex' != `numberOfArgs`function'Label`callIndex''
              #write <sector_`sectorID'_`cppOrder'.cpp> ","
            #EndIf
          #EndDo
          #write <sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
        #EndIf
      #EndIf

      multiply replace_(SecDecInternal`function'Call`callIndex'Unparsed,0 , SecDecInternalLabel`function'Call`callIndex'Arg,0 );
      .sort

    #EndDo

    #clearoptimize
*   }

    drop parseNext;
  #EndDo

  .sort
  drop unparsed;

* undefine largestLabel`...' after parsing
  #Do function = {`functionsToReplace'}
    #undefine largestLabel`function'
  #EndDo

* }

* Optimize the sign check (if applicable) and the final expression
* {
* reduce the labels SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheck`signCheckId' to one single label
  #If `contourDeformation'
    #Do signCheckId = 1, `numberOfRequiredSignChecks'
      Id SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheck`signCheckId' = SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheckGlobal ^ `signCheckId';
    #EndDo
  #EndIf

  Format float 20;
  Format C;
  Format O`optimizationLevel';
  #If `contourDeformation'
    Bracket SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheckGlobal;
  #EndIf
  .sort
  ExtraSymbols,array,SecDecInternalSecondAbbreviation;
  #optimize toOptimize

  #If `optimmaxvar_' > `numberOfSecondAbbreviations'
    #redefine numberOfSecondAbbreviations "`optimmaxvar_'"
  #EndIf

  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@##@SecDecInternalNewline@#// begin final dependency level#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
* }

* write code for the sign check
* {
  #If `contourDeformation'

    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#// contour deformation sign checks#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternalSignCheckExpression;#@SecDecInternalNewline@#"

    hide; nhide toOptimize, expr;
    .sort

    #Do signCheckId = 1, `numberOfRequiredSignChecks'

      Bracket SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheckGlobal;
      .sort
      Local expr = toOptimize[SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheckGlobal ^ `signCheckId'];
      .sort

      #write <sector_`sectorID'_`cppOrder'.cpp> "SecDecInternalSignCheckExpression = SecDecInternalImagPart(%%E);#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
      #write <sector_`sectorID'_`cppOrder'.cpp> "if (SecDecInternalSignCheckExpression > 0)"
      #write <sector_`sectorID'_`cppOrder'.cpp> "throw secdecutil::sign_check_error(#@SecDecInternalDblquote@##@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"

    #EndDo

    multiply replace_(SecDecInternalLabel`SecDecInternalContourDeformationPolynomial'CallSignCheckGlobal,0);
    .sort

    unhide;
    drop expr;
    .sort

    #write <sector_`sectorID'_`cppOrder'.cpp> "// end of contour deformation sign checks#@SecDecInternalNewline@#"

  #EndIf
* }

* write the integrand
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" toOptimize(#@no_split_expression@#)
  #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* return statement
  #write <sector_`sectorID'_`cppOrder'.cpp> "return tmp;#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables',sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`realParameters',sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`complexParameters',sector_`sectorID'_`cppOrder'.cpp)
  #If `contourDeformation'
    #call cppUndefine(`occurringDeformationParameters',sector_`sectorID'_`cppOrder'.cpp)
  #EndIf
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalDenominator#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalRealPart#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalImagPart#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalAbbreviations`shiftedOrderIndex'#@SecDecInternalNewline@#"

* Close the c++ function and namespaces
  #write <sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* Write the corresponding header "sector_`sectorID'_`cppOrder'.hpp".
  #write <sector_`sectorID'_`cppOrder'.hpp> "#ifndef `name'_codegen_sector_`sectorID'_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#define `name'_codegen_sector_`sectorID'_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#include \"`name'.hpp\"#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.hpp> "#include \"functions.hpp\"#@SecDecInternalNewline@#"
  #If `contourDeformation'
    #write <sector_`sectorID'_`cppOrder'.hpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #EndIf
  #write <sector_`sectorID'_`cppOrder'.hpp> "#define sector_`sectorID'_order_`cppOrder'_optimmaxvar_second `numberOfSecondAbbreviations'#@SecDecInternalNewline@#"
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
