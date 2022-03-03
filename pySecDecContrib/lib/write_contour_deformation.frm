* This file is included in "write_integrand.frm".
* It should be seen as a procedure that writes
* additional code for the contour deformation,
* namely the optimization of the parameters
* and the contour deformation polynomial.
* {

  #clearoptimize
  drop;
  Format normal;
  Format rational;
  Format 255;
  .sort

* define the argument for "replace_" that sets all regulators to zero
  #redefine nullifyRegulators ""
  #$counter = 1;
  #Do var = {`regulators'}
    #If `$counter' != 1
      #redefine nullifyRegulators "`nullifyRegulators' , "
    #EndIf
    #redefine nullifyRegulators "`nullifyRegulators' `var',0"
    #$counter = $counter + 1;
  #EndDo

* Define the argument for "replace_" that sets all absent integration variables to zero.
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

* Define the integration variables, the real parameters, the complex parameters,
* and the deformation parameters as c preprocessor variables.
* (The c function takes them packed into an array).
* "Format rational": Need the indices as integers.
  Format rational;
  Format 255;

* Optmimize the 2nd Symnanzik polynomial "F"

  #$counter = -1;
  #redefine deformedIntegrationVariables ""
  #Do IV = {`integrationVariables',}
    #If x`IV' != x
      #$counter = $counter + 1;
      #If `$counter' != 0
        #redefine deformedIntegrationVariables "`deformedIntegrationVariables',"
      #EndIf
      #redefine deformedIntegrationVariables "`deformedIntegrationVariables'SecDecInternalDeformed`IV'"
    #EndIf
  #EndDo

  Local expressionF = `SecDecInternalContourDeformationPolynomial'(`deformedIntegrationVariables',`regulators')*replace_(`nullifyRegulators')

* replace the calls to the deformed integration variables by dummy symbols
* {
  #Do IV = {`occurringIntegrationVariables',}
    #If x`IV' != x
      + SecDecInternalLabelSecDecInternalDeformed`IV' * SecDecInternalDeformed`IV'(`integrationVariables')
    #EndIf
  #EndDo
  ;
  argument `SecDecInternalContourDeformationPolynomial';
    #Do IV = {`occurringIntegrationVariables',}
      #If x`IV' != x
        Id SecDecInternalDeformed`IV' = SecDecInternalSecDecInternalDeformed`IV'Call;
      #EndIf
    #EndDo
* }
* remove calls to the deformation of absent integration variables
    #Do IV = {`absentIntegrationVariables',}
      #If x`IV' != x
        Id SecDecInternalDeformed`IV' = 0;
      #EndIf
    #EndDo
  endArgument;

* insert the deformation and `SecDecInternalContourDeformationPolynomial'
  #call insertDeformedIntegrationVariables
  argument SecDecInternalRealPart;
    #call insertDeformedIntegrationVariables
    #call insertDecomposed
  endArgument;
  #call insertDeformedIntegrationVariables
  #call insertDecomposed
  multiply replace_(
                        I,i_
                        #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
                          ,`nullifyAbsentIVs'
                        #EndIf
                   );
  .sort

* replace the calls to "SecDecInternalRealPart" by dummy symbols
* {
  #If `contourDeformation'

    #redefine function "SecDecInternalRealPart"
    #$labelCounter = 0;

    #Do i = 1,1

*     set dollar variable
      #$unmatched = 1;
      if ( $unmatched );
        #Do depth = 0, `insertionDepth'
          #call beginArgumentDepth(`depth')
            if ( $unmatched );
              if ( match(once `function'(SecDecInternalsDUMMY?$argLocal)) );
                $unmatched = 0;
                $arg = $argLocal;
                redefine i "0";
              endif;
            endif;
          #call endArgumentDepth(`depth')
        #EndDo
      endif;
      ModuleOption,sum,$arg;
      ModuleOption,minimum,$unmatched;
      ModuleOption,local,$argLocal;
      .sort:match;

*     The following "#if" evaluates to true only if there is still something to do.
      #If `i' == 0

        #$labelCounter = $labelCounter + 1;

        #Do depth = 0, `insertionDepth'
          Id `function'($arg) = SecDecInternal`function'Call`$labelCounter';
        #EndDo

        .sort

        Local expressionF = expressionF + SecDecInternalLabel`function'Call`$labelCounter' * (`$arg');

      #EndIf

    #EndDo

    #redefine largestLabel`function' "`$labelCounter'"

  #EndIf
* }

  Format O`optimizationLevel';
  AntiBracket `integrationVariables', `realParameters', `complexParameters', `deformationParameters';
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviation;

  #optimize expressionF

  Format float 20;
  Format C;
  Format 255;

* Write the function to optimize the contour deformation parameters
  #write <sector`sectorID'.info> "@order`shiftedOrderIndex'_contourDeformationPolynomialBody="
  #write <sector`sectorID'.info> "%O"

* c++ define the calls to "SecDecInternalRealPart"
  #redefine function "SecDecInternalRealPart"
  #Do labelID = 1, `largestLabel`function''
    Bracket SecDecInternalLabel`function'Call`labelID';
    .sort
    L realPart = expressionF[SecDecInternalLabel`function'Call`labelID'];
    .sort
    #write <sector`sectorID'.info> "SecDecInternal`function'Call`labelID' = SecDecInternalRealPart(%E);" realPart(#@FAIL@#)
    multiply replace_(SecDecInternalLabel`function'Call`labelID', 0);
    .sort
  #EndDo
  #undefine largestLabel`function'
  drop realPart;

* c++ define the defomed integration variables
  #Do IV = {`occurringIntegrationVariables',}
    #If x`IV' != x

      Bracket SecDecInternalLabelSecDecInternalDeformed`IV';
      .sort
      L deformedIV = expressionF[SecDecInternalLabelSecDecInternalDeformed`IV'];
      .sort
      #write <sector`sectorID'.info> "SecDecInternalSecDecInternalDeformed`IV'Call = %E;" deformedIV(#@FAIL@#)
      multiply replace_(SecDecInternalLabelSecDecInternalDeformed`IV', 0);
      .sort

    #EndIf
  #EndDo
  drop deformedIV;

* write `SecDecInternalContourDeformationPolynomial'
  #write <sector`sectorID'.info> "return(%E);" expressionF(#@FAIL@#)
  #write <sector`sectorID'.info> "@end"

* Delete "expressionF" since it is no longer needed.
  #clearoptimize
  drop expressionF;
  .sort


* The code above implements the contour deformation. Below, we implement an
* algorithm that finds reasonable `deformation_parameters`.


* The header can already completely be written here:
  #write <sector`sectorID'.info> "@order`shiftedOrderIndex'_optimizeDeformationParametersBody="

* The contour deformation is defined as:
* ``z_k({x_k}) = x_k - i * lambda_k * x_k * (1-x_k) * Re(dF_dx_k)``, where "dF_dx_k" denotes the derivative of `F` by `x_k`
* We compute ``1 / (  x_k * (1-x_k) * Re(dF_dx_k)  )`` for multiple ``x_k`` and - in the class "SectorContainerWithoutDeformation" -
* set `lambda_k` to the maximum in order to have the deformation smaller than one.

* Construct the ``  x_k * (1-x_k) * dF_dx_k  ``
* Note: We can take the real part of the full expression since ``x_k * (1-x_k)`` is real anyway.
  Local deformations =
  #Do idx = {`occurringIntegrationVariableIndices',}
    #If x`idx' != x
      + SecDecInternalLabelDeformation^(`idx'+1) *
        (
            `integrationVariable`idx'' * (1-`integrationVariable`idx'') *
            d`SecDecInternalContourDeformationPolynomial'd`idx'(`integrationVariables',`regulators')
        )
    #EndIf
  #EndDo
  ;

* Explicitly insert the derivatives of `SecDecInternalContourDeformationPolynomial'
  #call insertDecomposed

* Set the following variables to zero:
*   - all regulators
*   - absent integration variables (if any)
* translate sympy's imaginary unit to FORM's imaginary unit
  multiply replace_(
                       `nullifyRegulators'
                       #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
                         ,`nullifyAbsentIVs'
                       #EndIf
                       ,I,i_
                   );

* Define the integration variables, the real parameters, the complex parameters,
* the deformation parameters, "SecDecInternalRealPart", and "SecDecInternalAbs" as c preprocessor variables.
* (The c function takes them packed into an array).
* "Format rational": Need the indices as integers.
  Format rational;
  Format 255;

* optimize
  Bracket SecDecInternalLabelDeformation;
  Format O`optimizationLevel';
  Format 255;
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviation;
  #optimize deformations
  intohide deformations;
  Bracket SecDecInternalLabelDeformation;
  .sort

* write the optimization symbols
  Format float 20;
  Format C;
  Format 255;
  #write <sector`sectorID'.info> "%O"

* set the output lambdas to ``1 / Abs(    Re(  x_k * (1-x_k) * dF_dx_k  )    )``
  #$cppidx = -1;
  #Do idx = {`occurringIntegrationVariableIndices',}
    #If x`idx' != x
      #$cppidx = $cppidx + 1;
      L expr = deformations[SecDecInternalLabelDeformation^(`idx'+1)];
      .sort
      Format rational;
      Format C;
      Format 255;
      #write <sector`sectorID'.info> "SecDecInternalOutputDeformationParameters(`$cppidx',"
      Format float 20;
      Format C;
      Format 255;
      #write <sector`sectorID'.info> "1.0/SecDecInternalAbs(SecDecInternalRealPart(%E)));" expr(#@FAIL@#)
    #EndIf
  #EndDo

* Close the function and the namespace in the c++ file
  #write <sector`sectorID'.info> "@end"

* Delete the expressions "expr" and "deformations" since they are no longer needed.
  #clearoptimize
  unhide deformations;
  drop expr;
  drop deformations;
  .sort

* }
