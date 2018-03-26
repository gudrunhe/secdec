#-
Off statistics;

* Workaround:
* Do not use gzip since it sometimes fails for large expressions
* due to a bug in FORM.
On compress;

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

*Define a procedure that removes expressions known to be zero
#redefine knownZeroFunctions ""
#procedure nullify
  #Do depth = 0, `insertionDepth'
    #call beginArgumentDepth(`depth')
      #Do zeroFunction = {`knownZeroFunctions',}
        #If x`zeroFunction' != x
          Id `zeroFunction'(?args) = 0;
        #EndIf
      #EndDo
    #call endArgumentDepth(`depth')
  #EndDo
#endProcedure

* Define a simplification procedure
#procedure simplify
  #call nullify
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
      repeat Id SecDecInternalsDUMMY1? ^ (SecDecInternalsDUMMY2?!int_) = SecDecInternalPow(SecDecInternalsDUMMY1,SecDecInternalsDUMMY2);
      repeat Id SecDecInternalsDUMMY1? ^ (SecDecInternalsDUMMY2?neg_) = SecDecInternalDenominator(SecDecInternalsDUMMY1) ^ (-SecDecInternalsDUMMY2);
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

* Find functions that are equal to zero regardless of their arguments
#redefine zeros ""
#Do regulator = {`regulators'}
  #redefine zeros "`zeros',0"
#EndDo
#redefine groupOfFunctionsToConsider "decomposedPolynomialDerivatives,functions"
#If `contourDeformation'
  #redefine groupOfFunctionsToConsider "`groupOfFunctionsToConsider',deformedIntegrationVariableDerivativeFunctions"
#EndIf
#Do functionsToConsider = {`groupOfFunctionsToConsider'}
  #If `functionsToConsider' == decomposedPolynomialDerivatives
    #redefine insertProcedure "insertDecomposed"
  #ElseIf `functionsToConsider' == functions
    #redefine insertProcedure "insertOther"
  #Else
    #redefine insertProcedure "insertDeformedIntegrationVariables"
  #EndIf

  #$counter = 0;
  #Do function = {``functionsToConsider''}
    Local zeroCheck = `function'(`integrationVariables'`zeros');

    #call `insertProcedure'
    multiply replace_(I,i_);
    .sort

    #If termsin(zeroCheck)
      #$counter = $counter + 1;
      #If `$counter' == 1
        #redefine knownNonzeroFunctionsThisGroup "`function'"
      #Else
        #redefine knownNonzeroFunctionsThisGroup "`knownNonzeroFunctionsThisGroup',`function'"
      #EndIf
    #Else
      #redefine knownZeroFunctions "`knownZeroFunctions',`function'"
    #EndIf
  #EndDo
* ensure at least two functions (nasty edge cases for #Do loops otherwise)
  #If `$counter' == 0
    #redefine `functionsToConsider' "SecDecInternalfDUMMYNONEXISTENT1,SecDecInternalfDUMMYNONEXISTENT2"
  #ElseIf `$counter' == 1
    #redefine `functionsToConsider' "`knownNonzeroFunctionsThisGroup',SecDecInternalfDUMMYNONEXISTENT"
  #Else
    #redefine `functionsToConsider' "`knownNonzeroFunctionsThisGroup'"
  #EndIf
#EndDo
drop zeroCheck;

#call defineExpansion

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
Bracket `regulators';

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

* Find the calls to the contour deformation polynomials that need a sign check.
* The calls differ only in which Feynman paramters are set to zero or one. We
* investigate that by looking at the calls to "SecDecInternalCalI" (and its derivatives).
* {
  #If `contourDeformation'

    #$labelCounterF = 0;
    #$labelCounterU = 0;
    Local tmp = expression;
    Local signCheck = 0;
    .sort
    hide; nhide tmp, signCheck;

    #Do function = {`calIDerivatives'}

      #Do i = 1,1
*       set dollar variable
        skip signCheck;
        #$unmatched = 1;
        if ( $unmatched );
          if ( occurs(`function') );
            if ( match(once `function'(?SecDecInternalsDUMMY$args)) );
              redefine i "0";
              $wrappedArgs = SecDecInternalfDUMMYArgContainer($args);
            endif;
          endif;
        endif;
        ModuleOption,sum,$wrappedArgs;
        ModuleOption,local,$args;
        ModuleOption,local,$unmatched;
        .sort:match1;

*       The following "#if" evaluates to true only if there is still something to do.
        #If `i' == 0

          skip; nskip matcher;
          Local matcher = $wrappedArgs;
          Id once SecDecInternalfDUMMYArgContainer(?arguments$args) = 0;
          .sort:match2;
          drop matcher;

*         If all args are zero or one, the expression to be sign checked evaluates to zero
*         --> can omit that check because it always passes.
          #redefine needThisCheck "0"
          #Do arg = {`$args',}
            #If x`arg' != x
              #If (`arg' != 0) && (`arg' != 1)
                #redefine needThisCheck "1"
              #EndIf
            #EndIf
          #EndDo

          #If `needThisCheck'

            #$labelCounterF = $labelCounterF + 1;

*           The sign check to be performed for the contour deformation polynomial (F) is:
*           "SecDecInternalImagPart( SecDecInternalContourDeformationPolynomial(<deformed integration variables>) - SecDecInternalContourDeformationPolynomial(<undeformed integration variables>) ) <= 0"
            skip;
            Local signCheck = signCheck +
                SecDecInternalLabelContourDeformationPolynomialCallSignCheck`$labelCounterF' *
                (
                  SecDecInternalfDUMMYdeformedContourDeformationPolynomial($args) - `SecDecInternalContourDeformationPolynomial'($args)
                );
            .sort

*           The sign check to be performed for the positive polynomials (e.g. U) is:
*           "SecDecInternalRealPart( SecDecInternalPositivePolynomial(<deformed integration variables>) ) >= 0"
            #Do positivePolynomial = {`positivePolynomials',}
              #If x`positivePolynomial' != x
                #$labelCounterU = $labelCounterU + 1;
                skip;
                Local signCheck = signCheck +
                    SecDecInternalLabelUCallSignCheck`$labelCounterU' * SecDecInternalfDUMMY`positivePolynomial'($args);
                .sort
              #EndIf
            #EndDo

*           Define a variables that only keeps the integration variables in "$args"
*           {
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
*           }

            #Do IV = {`integrationVariables'}
              Id SecDecInternalfDUMMYdeformedContourDeformationPolynomial(?frontArgs,`IV',?backArgs) =
                 SecDecInternalfDUMMYdeformedContourDeformationPolynomial(?frontArgs,SecDecInternalDeformed`IV'(`argsWithoutRegulators'),?backArgs);

              #Do positivePolynomial = {`positivePolynomials',}
                #If x`positivePolynomial' != x
                  Id SecDecInternalfDUMMY`positivePolynomial'(?frontArgs,`IV',?backArgs) =
                     SecDecInternalfDUMMY`positivePolynomial'(?frontArgs,SecDecInternalDeformed`IV'(`argsWithoutRegulators'),?backArgs);
                #EndIf
              #EndDo
            #EndDo

            Id SecDecInternalfDUMMYdeformedContourDeformationPolynomial(?args) = `SecDecInternalContourDeformationPolynomial'(?args);

            #Do positivePolynomial = {`positivePolynomials',}
              #If x`positivePolynomial' != x
                Id SecDecInternalfDUMMY`positivePolynomial'(?args) = `positivePolynomial'(?args);
              #EndIf
            #EndDo
          #EndIf

          #Do innerFunction = {`calIDerivatives'}
            Id `innerFunction'(`$args') = 0;
          #EndDo
          .sort

        #EndIf
      #EndDo

      #redefine numberOfRequiredFSignChecks "`$labelCounterF'"
      #redefine numberOfRequiredUSignChecks "`$labelCounterU'"

    #EndDo

    unhide;
    drop tmp;

  #EndIf
* }

  Local toOptimize = 0;

* insert calI
  #call insertCalI
  #call nullify
  .sort

* Replace calls to the deformation's Jacobian determinant "SecDecInternalContourdefJacobian" by symbols.
* {

  #If `contourDeformation'

*   Make sure all labels exist.
    #Do function = {`contourdefJacobianFunctions',}
      #If x`function' != x
        #redefine largestLabel`function' "0"
      #EndIf
    #EndDo

*   No need for "#Do depth = 0, `insertionDepth'" because the Jacobian determinant should only occur at top level.

*   Since we need intermediate ".sort" instructions, we cannot use the
*   "repeat" environment.
*   The following construction is suggested in the FORM documentation.

    skip toOptimize;
    #Do i = 1,1

*     set dollar variable
      #$unmatched = 1;
      if ( $unmatched );
        #Do depth = 0, `insertionDepth'
          #call beginArgumentDepth(`depth')
            if ( $unmatched );
              if ( match(once SecDecInternalfDUMMY?{`contourdefJacobianFunctions'}$function(?SecDecInternalsDUMMY$args)) );
                $unmatched = 0;
                $wrappedFunctionWithArgs = SecDecInternalfDUMMYwrappedFunctionWithArgs($function, $args);
                redefine i "0";
              endif;
            endif;
          #call endArgumentDepth(`depth')
        #EndDo
      endif;
      ModuleOption,sum,$wrappedFunctionWithArgs;
      ModuleOption,minimum,$unmatched;
      ModuleOption,local,$function;
      ModuleOption,local,$args;
      .sort:match1;

*     The following "#if" evaluates to true only if there is still something to do.
      #If `i' == 0

        skip; nskip matcher;
        Local matcher = $wrappedFunctionWithArgs;
        Id once SecDecInternalfDUMMYwrappedFunctionWithArgs(SecDecInternalfDUMMY?$function, ?arguments$args) = 0;
        .sort:match2;
        drop matcher;

        #IfDef `largestLabel`$function''
          #$labelCounter = `largestLabel`$function'' + 1;
        #Else
          #$labelCounter = 1;
        #EndIf

        skip;
        Local toOptimize = toOptimize +
            SecDecInternalLabel`$function'Call`$labelCounter'Arg * `$function'($args);
        .sort
        skip toOptimize;

        Id `$function'(`$args') = SecDecInternal`$function'Call`$labelCounter';

*       no arguments after "insertContourdefJacobianDerivatives"
        #redefine numberOfArgs`$function'Label`$labelCounter' "0"

        #redefine largestLabel`$function' "`$labelCounter'"

      #EndIf
    #EndDo

*   insert the contour deformation Jacobian and its derivatives
    #call insertContourdefJacobianDerivatives


  #EndIf

* }

* Explicitly insert the functions defined in python.
* {

  #If `contourDeformation'
    #redefine insertProceduresToConsider "insertDeformedIntegrationVariables,insertOther,insertDeformedIntegrationVariables,insertDecomposed"
  #Else
    #redefine insertProceduresToConsider "insertOther,insertDecomposed"
  #EndIf

  #Do insertProcedure = {`insertProceduresToConsider'}

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
        repeat Id SecDecInternalDenominator(SecDecInternalsDUMMY?number_) = 1/SecDecInternalsDUMMY;
        repeat Id SecDecInternalfDUMMY?(?SecDecInternalsDUMMY) * SecDecInternalDenominator(SecDecInternalfDUMMY?(?SecDecInternalsDUMMY)) = 1;

*       Remove vanishing deformed variables
        #If `contourDeformation'
          #call removeVanishingDeformedIntegrationVariableDerivatives
        #EndIf

      #call endArgumentDepth(`depth')

    #EndDo

*   Replace calls to the known functions (and appearing derivatives) by symbols.
*   {

    #If `insertProcedure' == insertOther
      #redefine functionsToReplace "`functions'"
    #ElseIf `insertProcedure' == insertDeformedIntegrationVariables
      #redefine functionsToReplace "`deformedIntegrationVariableDerivativeFunctions'"
    #Else
      #redefine functionsToReplace "`decomposedPolynomialDerivatives'"
    #EndIf

*   Make sure all labels exist.
    #Do function = {`functionsToReplace'}
      #IfNDef `largestLabel`function''
        #redefine largestLabel`function' "0"
      #EndIf
    #EndDo

*   Since we need intermediate ".sort" instructions, we cannot use the
*   "repeat" environment.
*   The following construction is suggested in the FORM documentation.

    #Do i = 1,1
      Bracket `functionsToReplace';
      .sort:bracket;
*     set dollar variables
      #$unmatched = 1;
      if ( $unmatched );
        #Do depth = 0, `insertionDepth'
          #call beginArgumentDepth(`depth')
            if ( $unmatched );
              if ( match(once SecDecInternalfDUMMY?{`functionsToReplace'}$function(?SecDecInternalsDUMMY$args)) );
                $unmatched = 0;
                $wrappedFunctionWithArgs = SecDecInternalfDUMMYwrappedFunctionWithArgs($function, $args);
                redefine i "0";
              endif;
            endif;
          #call endArgumentDepth(`depth')
        #EndDo
      endif;
      ModuleOption,sum,$wrappedFunctionWithArgs;
      ModuleOption,minimum,$unmatched;
      ModuleOption,local,$function;
      ModuleOption,local,$args;
      .sort:match;

*     The following "#if" evaluates to true only if there is still something to do.
      #If `i' == 0

        skip; nskip toUnwrap;
        Local toUnwrap = $wrappedFunctionWithArgs;
        Id,once SecDecInternalfDUMMYwrappedFunctionWithArgs(SecDecInternalfDUMMY?$function, ?args$args) = SecDecInternalfDUMMY(?args);
        #call `insertProcedure'
        #call nullify
        .sort:unwrap;
        drop toUnwrap;

        #If termsin(toUnwrap) == 0

          if ( occurs(`$function') );
            #Do depth = 0, `insertionDepth'
              #call beginArgumentDepth(`depth')
                if ( occurs(`$function') ) Id `$function'($args) = 0;
              #call endArgumentDepth(`depth')
            #EndDo
          endif;

        #ElseIf termsin(toUnwrap) == 1

          if ( occurs(`$function') );
            #Do depth = 0, `insertionDepth'
              #call beginArgumentDepth(`depth')
                if ( occurs(`$function') ) Id `$function'($args) = SecDecInternalfDUMMYContainer(SecDecInternalsDUMMY`$function',$args);
              #call endArgumentDepth(`depth')
            #EndDo
          endif;

        #Else

          #$labelCounter = `largestLabel`$function'' + 1;
          #redefine largestLabel`$function' "`$labelCounter'"

*         no arguments after "insertDeformedIntegrationVariables"
          #redefine numberOfArgs`$function'Label`$labelCounter' "0"

          Local toOptimize = toOptimize + SecDecInternalLabel`$function'Call`$labelCounter'Arg * SecDecInternalfDUMMYContainer(SecDecInternalsDUMMY`$function',`$args');

          if ( occurs(`$function') );
            #Do replaceDepth = 0, `insertionDepth'
              #call beginArgumentDepth(`replaceDepth')
                if ( occurs(`$function') ) Id `$function'($args) = SecDecInternal`$function'Call`$labelCounter';
              #call endArgumentDepth(`replaceDepth')
            #EndDo
          endif;

        #EndIf

      #EndIf

    #EndDo

    #Do function = {`functionsToReplace'}
      #Do depth = 0, `insertionDepth'
        #call beginArgumentDepth(`depth')
          repeat Id SecDecInternalfDUMMYContainer(SecDecInternalsDUMMY`function',?arguments) = `function'(?arguments);
        #call endArgumentDepth(`depth')
      #EndDo
    #EndDo
    .sort:backsubs;

*   }

    #Do depth = 0, `insertionDepth'
      #call beginArgumentDepth(`depth')
        #call `insertProcedure'
        #call nullify
      #call endArgumentDepth(`depth')
    #EndDo
    .sort

*   some simplifications
    #Do depth = 0, `insertionDepth'
      #call beginArgumentDepth(`depth')
        repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 0) = 1;
        repeat Id SecDecInternalPow(SecDecInternalsDUMMYbase?, 1) = SecDecInternalsDUMMYbase;
        repeat Id SecDecInternalPow(0, SecDecInternalsDUMMYexponent?) = 0;
        Denominators SecDecInternalDenominator;
        factarg,(-1),SecDecInternalDenominator;
        chainout SecDecInternalDenominator;
        repeat Id log(1) = 0;
        repeat Id SecDecInternalsDUMMY1? ^ (SecDecInternalsDUMMY2?!int_) = SecDecInternalPow(SecDecInternalsDUMMY1,SecDecInternalsDUMMY2);
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

  #EndDo

* }

* translate sympy's imaginary unit to FORM's imaginary unit
  multiply replace_(I,i_);
  .sort

* Find and count the occurring integration variables.
* {
* "Format rational" because we need the dollar variables as integers
  Format rational;
  hide; nhide toOptimize, expression;
  #call simplify

  #redefine occurringIntegrationVariables ""
  #redefine occurringIntegrationVariableIndices ""
  #redefine absentIntegrationVariables ""
  #redefine absentIntegrationVariableIndices ""
  #$counterOccur = 0;
  #$counterAbsent = 0;
  #$currentIVIndex = -1;

  #Do IV = {`integrationVariables'}
    #redefine `IV'Occurs "0"
    if ( occurs(`IV') ) redefine `IV'Occurs "1";
  #EndDo
  .sort

  #Do IV = {`integrationVariables'}

    #$currentIVIndex = $currentIVIndex + 1;

    #If ``IV'Occurs'
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

* Explicitly insert calls with only one term since these do not lead
* to an undesired expansion of expressions.
* {

  #If `contourDeformation'
* `decomposedPolynomialDerivatives' appears twice to account for "decomposedPoly(deformedVariable)"
  #redefine functionsToReplace "`decomposedPolynomialDerivatives',`deformedIntegrationVariableDerivativeFunctions',`decomposedPolynomialDerivatives',`functions'"
  #Else
  #redefine functionsToReplace "`decomposedPolynomialDerivatives',`functions'"
  #Endif

  #Do function = {`functionsToReplace'}
    #Do callIndex = 1, `largestLabel`function''

      Bracket SecDecInternalLabel`function'Call`callIndex'Arg;
      .sort

      skip; nskip thisCall;
      L thisCall = toOptimize[SecDecInternalLabel`function'Call`callIndex'Arg];
      #call simplify
      .sort
      drop thisCall;

      #If termsin(thisCall) <= 1
        multiply replace_(SecDecInternal`function'Call`callIndex',thisCall);
      #EndIf

    #EndDo
  #EndDo
  .sort

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

* Since we need intermediate ".sort" instructions, we cannot use the
* "repeat" environment.
* The following construction is suggested in the FORM documentation.

  #Do i = 1,1

*   set dollar variables
    #$unmatched = 1;
    if ( $unmatched );
      #Do depth = 0, `insertionDepth'
        #call beginArgumentDepth(`depth')
          if ( $unmatched );
            if ( match(once SecDecInternalfDUMMY?{`functionsToReplace'}$function(?SecDecInternalsDUMMY$args)) );
              $unmatched = 0;
              $wrappedfunctionWithArgs = SecDecInternalfDUMMYwrappedFunctionWithArgs($function, $args);
              redefine i "0";
            endif;
          endif;
        #call endArgumentDepth(`depth')
      #EndDo
    endif;
    ModuleOption,sum,$wrappedfunctionWithArgs;
    ModuleOption,minimum,$unmatched;
    ModuleOption,local,$function;
    ModuleOption,local,$args;
    .sort:match1;

*   The following "#if" evaluates to true only if there is still something to do.
    #If `i' == 0

      skip; nskip matcher;
      Local matcher = $wrappedfunctionWithArgs;
      Id once SecDecInternalfDUMMYwrappedFunctionWithArgs(SecDecInternalfDUMMY?$function, ?arguments$args) = 0;
      .sort:match2;
      drop matcher;

      #$labelCounter = `largestLabel`$function'' + 1;

      if ( occurs(`$function') );
        #Do replaceDepth = 0, `insertionDepth'
          #call beginArgumentDepth(`replaceDepth')
            Id `$function'($args) = SecDecInternal`$function'Call`$labelCounter';
          #call endArgumentDepth(`replaceDepth')
        #EndDo
      endif;
      .sort

      skip; nskip arguments;
      L arguments = SecDecInternalfDUMMYarguments($args);

      repeat Id SecDecInternalfDUMMYarguments(SecDecInternalsDUMMY?, ?otherArgs) =
          SecDecInternalLabel`$function'Call`$labelCounter'Arg * (SecDecInternalsDUMMY + SecDecInternalfDUMMYarguments(?otherArgs));

*     Define `$argCounter' by loking at the term with the empty function "SecDecInternalfDUMMYarguments"
      Id SecDecInternalfDUMMYarguments * SecDecInternalLabel`$function'Call`$labelCounter'Arg ^ SecDecInternalsDUMMYexponent?$argCounter = 0;

      .sort
      drop arguments;

*     Add all arguments to top level polynomial for simultaneous optimization.
      Local toOptimize = toOptimize + arguments;

      #redefine numberOfArgs`$function'Label`$labelCounter' "`$argCounter'"
      #redefine largestLabel`$function' "`$labelCounter'"

    #EndIf

  #EndDo

  drop expression;
  #If `contourDeformation'
    drop signCheck;
  #EndIf

  Local toOptimize = toOptimize + expression
  #If `contourDeformation'
     + signCheck
  #EndIf
  ;
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
  #write <sector_`sectorID'_`cppOrder'.cpp> "#ifdef SECDEC_WITH_CUDA#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) (complex_t{x}).real()#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalImagPart(x) (complex_t{x}).imag()#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#else#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) std::real(x)#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalImagPart(x) std::imag(x)#@SecDecInternalNewline@#"
  #write <sector_`sectorID'_`cppOrder'.cpp> "#endif#@SecDecInternalNewline@#"

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
  AutoDeclare CFunctions SecDecInternalAbbreviations`shiftedOrderIndex';
  ExtraSymbols,array,SecDecInternalAbbreviations`shiftedOrderIndex';
  #$toOptimize = toOptimize;
  #clearoptimize
  Format O0;
  Format normal;
  Local toOptimize = `$toOptimize';

* Define the function calls replaced by symbols.
* The difficulty here is that the function arguments may
* refer to other function calls. --> We must order the
* definitions in the c++ file properly.
* {

* Keep track of function calls that are not written to
* the c++ file yet.

  #redefine functionsToReplace "`functions',`decomposedPolynomialDerivatives',log,SecDecInternalPow,SecDecInternalDenominator"
  #If `contourDeformation'
    #redefine functionsToReplace "SecDecInternalRealPart,`contourdefJacobianFunctions',`deformedIntegrationVariableDerivativeFunctions',`functionsToReplace'"
  #EndIf

  L unparsed = 0;
  #Do function = {`functionsToReplace'}
    #Do callIndex = 1, `largestLabel`function''
      .sort
      skip;
      L unparsed = unparsed + SecDecInternal`function'Call`callIndex'Unparsed * `function'(`callIndex');
    #EndDo
  #EndDo
  .sort
  skip unparsed;
  #Do function = {`functionsToReplace'}
    #Do callIndex = 1, `largestLabel`function''
      Id SecDecInternal`function'Call`callIndex' = SecDecInternal`function'Call`callIndex' * SecDecInternalsDUMMYhaveUnparsedDependencies;
    #EndDo
  #EndDo
  .sort

  #$globalLabelCounter = 1;
  #$globalIDCounter = 1;
  #redefine numberOfSecondAbbreviations "0"

  #Do i = 1,1
    Format normal;
    skip;
    .sort
    skip;
    #redefine parseNextGlobalIDMin "`$globalIDCounter'"
    #redefine parseNextGlobalIDMax "0"
    Local parseNext = 0;
    .sort

    #Do unparsedTerm = unparsed

      skip; nskip thisTerm;
      Local thisTerm = `unparsedTerm';
      Id SecDecInternalfDUMMY?$function(SecDecInternalsDUMMY?$callIndex) = 0;
      .sort:getTerm;
      drop thisTerm;

      skip; nskip toOptimize;
      Bracket SecDecInternalLabel`$function'Call`$callIndex'Arg;
      .sort

*     We can write the call under consideration to the c++ file only
*     if all dependent calls are already written.
      skip; nskip unparsed, expr, arg1,...,arg`numberOfArgs`$function'Label`$callIndex'';
      #If `numberOfArgs`$function'Label`$callIndex'' == 0
        Local expr = toOptimize[SecDecInternalLabel`$function'Call`$callIndex'Arg];
      #Else
        #Do argIndex = 1, `numberOfArgs`$function'Label`$callIndex''
           Local arg`argIndex' = toOptimize[SecDecInternalLabel`$function'Call`$callIndex'Arg ^ `argIndex'];
        #EndDo
      #EndIf
      #redefine dependenciesDone "1"
      if ( occurs(SecDecInternalsDUMMYhaveUnparsedDependencies) );
        redefine dependenciesDone "0";
        redefine i "0";
      endif;
      .sort

      #If `dependenciesDone'

*       Define the variable only if it is used.
        #redefine needThisCall "0"
        skip unparsed;
        if ( occurs(SecDecInternal`$function'Call`$callIndex') ) redefine needThisCall "1";
        .sort

*       Prepare for another optimization
*       {
        #If `needThisCall'

          #redefine parseNextGlobalIDMax "`$globalIDCounter'"
          #redefine parseNextFunction`$globalIDCounter' "`$function'"
          #redefine parseNextCallIndex`$globalIDCounter' "`$callIndex'"
          #$globalIDCounter = $globalIDCounter + 1;
          #redefine SecDecInternalLabel`$function'Call`$callIndex'GlobalIDMin "`$globalLabelCounter'"
          Format normal;
          skip; nskip parseNext;

          #If `numberOfArgs`$function'Label`$callIndex'' == 0
            Local parseNext = parseNext
              + expr * SecDecInternalLabelSecDecInternalGeneral ^ `$globalLabelCounter';
            #$globalLabelCounter = $globalLabelCounter + 1;
          #Else
            Local parseNext = parseNext
            #Do argIndex = 1, `numberOfArgs`$function'Label`$callIndex''
              + arg`argIndex' * SecDecInternalLabelSecDecInternalGeneral ^ (`$globalLabelCounter' + `argIndex' - 1)
            #EndDo
            ;
            #$globalLabelCounter = $globalLabelCounter + `numberOfArgs`$function'Label`$callIndex'';
          #EndIf

        #Else

          multiply replace_(SecDecInternal`$function'Call`$callIndex'Unparsed,0 , SecDecInternalLabel`$function'Call`$callIndex'Arg,0);

        #EndIf
*       }

      #EndIf

      .sort
      drop expr, arg1,...,arg`numberOfArgs`$function'Label`$callIndex'';

    #EndDo

*   optimize and write to c file
*   {
    skip unparsed;
    Format float 20;
    Format C;
    Format O`optimizationLevel';
    Bracket SecDecInternalLabelSecDecInternalGeneral;
    .sort
    ExtraSymbols,array,SecDecInternalSecondAbbreviation;
    skip unparsed;
    #optimize parseNext
    intohide parseNext;
    Bracket SecDecInternalLabelSecDecInternalGeneral;
    .sort

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
        skip unparsed;
        #write <sector_`sectorID'_`cppOrder'.cpp> "auto SecDecInternal`function'Call`callIndex' = "
        #If `numberOfArgs`function'Label`callIndex'' == 0
          Local expr = parseNext[SecDecInternalLabelSecDecInternalGeneral ^ `SecDecInternalLabel`function'Call`callIndex'GlobalIDMin'];
          .sort
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%e#@SecDecInternalNewline@#"  expr(#@no_split_expression@#)
          drop expr;
        #ElseIf `function' == SecDecInternalPow
          Local base = parseNext[SecDecInternalLabelSecDecInternalGeneral ^ (`SecDecInternalLabel`function'Call`callIndex'GlobalIDMin')];
          Local exponent = SecDecInternalfDUMMY(parseNext[SecDecInternalLabelSecDecInternalGeneral ^ (`SecDecInternalLabel`function'Call`callIndex'GlobalIDMin'+1)]);
          .sort
          #write <sector_`sectorID'_`cppOrder'.cpp> "pow(%%E," base(#@no_split_expression@#)
          drop base;
          skip; nskip exponent;
          #redefine exponentIsInteger "0"
          if ( match(once SecDecInternalfDUMMY(SecDecInternalsDUMMY?int_)) ) redefine exponentIsInteger "1";
          Id SecDecInternalfDUMMY(SecDecInternalsDUMMY?) = SecDecInternalsDUMMY;
          .sort
          #If `exponentIsInteger'
            Format rational;
          #EndIf
          #write <sector_`sectorID'_`cppOrder'.cpp> "%%E);#@SecDecInternalNewline@#" exponent(#@no_split_expression@#)
          Format C;
          Format float 20;
          drop exponent;
        #Else
          #write <sector_`sectorID'_`cppOrder'.cpp> "`function'("
          #Do argIndex = 1, `numberOfArgs`function'Label`callIndex''
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

        Id SecDecInternal`function'Call`callIndex' = SecDecInternal`function'Call`callIndex' / SecDecInternalsDUMMYhaveUnparsedDependencies;
        multiply replace_(SecDecInternal`function'Call`callIndex'Unparsed,0 , SecDecInternalLabel`function'Call`callIndex'Arg,0 );
        .sort

      #EndIf

    #EndDo

    #clearoptimize
*   }

    unhide parseNext;
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
* reduce the labels SecDecInternalLabel`ContourDeformationPolynomial/PositivePolynomial'CallSignCheck`signCheckId' to one two labels
  #If `contourDeformation'
    #Do signCheckId = 1, `numberOfRequiredFSignChecks'
      Id SecDecInternalLabelContourDeformationPolynomialCallSignCheck`signCheckId' = SecDecInternalLabelContourDeformationPolynomialCallSignCheckGlobal ^ `signCheckId';
    #EndDo
    #Do signCheckId = 1, `numberOfRequiredUSignChecks'
      Id SecDecInternalLabelUCallSignCheck`signCheckId' = SecDecInternalLabelUCallSignCheckGlobal ^ `signCheckId';
    #EndDo
  #EndIf

  Format float 20;
  Format C;
  Format O`optimizationLevel';
  #If `contourDeformation'
    Bracket SecDecInternalLabelContourDeformationPolynomialCallSignCheckGlobal, SecDecInternalLabelUCallSignCheckGlobal;
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

* write code for the sign checks
* {
  #If `contourDeformation'

    #write <sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#// contour deformation sign checks#@SecDecInternalNewline@#"
    #write <sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternalSignCheckExpression;#@SecDecInternalNewline@#"

    hide; nhide toOptimize, expr;
    Bracket SecDecInternalLabelContourDeformationPolynomialCallSignCheckGlobal;
    .sort

    #Do signCheckId = 1, `numberOfRequiredFSignChecks'

      keep brackets;
      skip; nskip expr;
      Local expr = toOptimize[SecDecInternalLabelContourDeformationPolynomialCallSignCheckGlobal ^ `signCheckId'];
      .sort

      #If termsin(expr) > 0
        #write <sector_`sectorID'_`cppOrder'.cpp> "SecDecInternalSignCheckExpression = SecDecInternalImagPart(%%E);#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
        #write <sector_`sectorID'_`cppOrder'.cpp> "#ifdef SECDEC_WITH_CUDA#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  if (SecDecInternalSignCheckExpression > 0) {"
        #write <sector_`sectorID'_`cppOrder'.cpp> "    printf(#@SecDecInternalDblquote@#Sign check `signCheckId' (contour deformation polynomial) failed.#@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "    return std::nan(#@SecDecInternalDblquote@##@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  }#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#else#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  if (SecDecInternalSignCheckExpression > 0)"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  throw secdecutil::sign_check_error(#@SecDecInternalDblquote@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  , #@SecDecInternalEscapedDblquote@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  contour deformation polynomial"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  #@SecDecInternalEscapedDblquote@#, check id #@SecDecInternalEscapedDblquote@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  `signCheckId'"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  #@SecDecInternalEscapedDblquote@#,#@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#endif#@SecDecInternalNewline@#"
      #EndIf

    #EndDo

    multiply replace_(SecDecInternalLabelContourDeformationPolynomialCallSignCheckGlobal,0);
    Bracket SecDecInternalLabelUCallSignCheckGlobal;
    .sort

    #Do signCheckId = 1, `numberOfRequiredUSignChecks'

      keep brackets;
      skip; nskip expr;
      Local expr = toOptimize[SecDecInternalLabelUCallSignCheckGlobal ^ `signCheckId'];
      .sort

      #If termsin(expr) > 0
        #write <sector_`sectorID'_`cppOrder'.cpp> "SecDecInternalSignCheckExpression = SecDecInternalRealPart(%%E);#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
        #write <sector_`sectorID'_`cppOrder'.cpp> "#ifdef SECDEC_WITH_CUDA#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  if (SecDecInternalSignCheckExpression < 0) {"
        #write <sector_`sectorID'_`cppOrder'.cpp> "    printf(#@SecDecInternalDblquote@#Sign check `signCheckId' (positive polynomial) failed.#@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "    return std::nan(#@SecDecInternalDblquote@##@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  }#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#else#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  if (SecDecInternalSignCheckExpression < 0)"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  throw secdecutil::sign_check_error(#@SecDecInternalDblquote@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  , #@SecDecInternalEscapedDblquote@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  positive polynomial"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  #@SecDecInternalEscapedDblquote@#, check id #@SecDecInternalEscapedDblquote@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  `signCheckId'"
        #write <sector_`sectorID'_`cppOrder'.cpp> "  #@SecDecInternalEscapedDblquote@#,#@SecDecInternalDblquote@#);#@SecDecInternalNewline@#"
        #write <sector_`sectorID'_`cppOrder'.cpp> "#endif#@SecDecInternalNewline@#"
      #EndIf

    #EndDo

    multiply replace_(SecDecInternalLabelUCallSignCheckGlobal,0);
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
* We still need a file that collects the whole sector.

* clear last step
.store

* "Format rational": Need the indices as integers.
Format rational;

* include series class
#write <sector_`sectorID'.cpp> "#include <secdecutil/series.hpp>#@SecDecInternalNewline@#"

#write <sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"

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
  #write <sector_`sectorID'.cpp> "#define sector_`sectorID'_order_`cppOrder'_numIV `numOccurringIVOrder`shiftedOrderIndex''#@SecDecInternalNewline@#"

* include the headers for all orders in this sector
  #write <sector_`sectorID'.cpp> "#include \"sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"

* include contour deformation and optimize deformation_parameter headers (if needed)
  #If `contourDeformation'
    #write <sector_`sectorID'.cpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
    #write <sector_`sectorID'.cpp> "#include \"optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #EndIf

  #write <sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"

#EndDo

* open c++ namespace
#write <sector_`sectorID'.cpp> "namespace `name'#@SecDecInternalNewline@#"
#write <sector_`sectorID'.cpp> "{#@SecDecInternalNewline@#"

* write the getter function of the container
#write <sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"
#write <sector_`sectorID'.cpp> "`integrandContainerType' get_integrand_of_sector_`sectorID'()#@SecDecInternalNewline@#"
#write <sector_`sectorID'.cpp> "{#@SecDecInternalNewline@#"
#write <sector_`sectorID'.cpp> "return `integrandContainerInitializer';#@SecDecInternalNewline@#"
#write <sector_`sectorID'.cpp> "};#@SecDecInternalNewline@#"
#write <sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"

* close c++ namespace
#write <sector_`sectorID'.cpp> "};#@SecDecInternalNewline@#"

#write <sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"

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
  #write <sector_`sectorID'.cpp> "#undef sector_`sectorID'_order_`cppOrder'_numIV#@SecDecInternalNewline@#"
#EndDo
*}

.end
