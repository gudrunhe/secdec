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
  .sort

* We are writing the c++ file "contour_deformation_sector_`sectorID'_`cppOrder'.cpp"
* and the corresponding header "contour_deformation_sector_`sectorID'_`cppOrder'.hpp".
* The header can already completely be written here:
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#ifndef `name'_codegen_contour_deformation_sector_`sectorID'_order_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#define `name'_codegen_contour_deformation_sector_`sectorID'_order_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include \"`name'.hpp\"#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include \"functions.hpp\"#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "{#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "  secdecutil::SectorContainerWithDeformation<real_t, complex_t>::DeformedIntegrandFunction sector_`sectorID'_order_`cppOrder'_contour_deformation_polynomial;#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "};#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#endif#@SecDecInternalNewline@#"

* Write includes and open the namespace in the corresponding .cpp file
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#include \"contour_deformation_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_complex_math.h>#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_linalg.h>#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "{#@SecDecInternalNewline@#"

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
  #call cppDefine(`occurringIntegrationVariables',integration_variables,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`realParameters',real_parameters,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`complexParameters',complex_parameters,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`occurringDeformationParameters',deformation_parameters,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)

* Since FORM does not use "abbreviation[0]", we can use it as temporary variable.
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

* c++-define "SecDecInternalRealPart"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) std::real(x)#@SecDecInternalNewline@#"

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

  L expressionF = `SecDecInternalContourDeformationPolynomial'(`deformedIntegrationVariables',`regulators')*replace_(`nullifyRegulators') + SecDecInternalsDUMMYExpressionFTail;
  .sort

* replace the calls to the deformed integration variables by dummy symbols
* {
  #Do IV = {`occurringIntegrationVariables',}
    #If x`IV' != x

      argument `SecDecInternalContourDeformationPolynomial';
        Id SecDecInternalDeformed`IV' = SecDecInternalSecDecInternalDeformed`IV'Call;
      endArgument;
      Id SecDecInternalsDUMMYExpressionFTail = SecDecInternalsDUMMYExpressionFTail +
        SecDecInternalLabelSecDecInternalDeformed`IV' * SecDecInternalDeformed`IV'(`integrationVariables');

    #EndIf
  #EndDo
* }

* insert the deformation and `SecDecInternalContourDeformationPolynomial'
  #call insertDeformedIntegrationVariables
  argument SecDecInternalRealPart;
    #call insertDeformedIntegrationVariables
    #call insertOther
  endArgument;
  #call insertDeformedIntegrationVariables
  #call insertOther
  multiply replace_(I,i_);
  .sort

* replace the calls to "SecDecInternalRealPart" by dummy symbols
* {
  #If `contourDeformation'

    #redefine function "SecDecInternalRealPart"
    #$labelCounter = 0;

    #Do i = 1,1
*     set dollar variable
      #Do depth = 0, `insertionDepth'
        if ( match(`function'(?SecDecInternalsDUMMY$arg)) ) redefine i "0";
      #EndDo
      .sort

*     The following "#if" evaluates to true only if there is still something to do.
      #If `i' == 0

        #$labelCounter = $labelCounter + 1;

        #Do depth = 0, `insertionDepth'
          Id `function'($arg) = SecDecInternal`function'Call`$labelCounter';
        #EndDo

        Id SecDecInternalsDUMMYExpressionFTail = SecDecInternalsDUMMYExpressionFTail +
            SecDecInternalLabel`function'Call`$labelCounter' * (`$arg');

      #EndIf
      .sort
    #EndDo

    #redefine largestLabel`function' "`$labelCounter'"

  #EndIf
* }

  Id SecDecInternalsDUMMYExpressionFTail = 0;

* remove calls to the deformation of absent integration variables
  #Do IV = {`absentIntegrationVariables',}
    #If x`IV' != x
      Id SecDecInternalDeformed`IV' = 0;
    #EndIf
  #EndDo


  #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
    multiply replace_(`nullifyAbsentIVs');
    .sort
  #EndIf

  Format O`optimizationLevel';
  AntiBracket `integrationVariables', `realParameters', `complexParameters', `deformationParameters';
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviation;

  #optimize expressionF

  Format float 20;
  Format C;

* Write the function to optimize the contour deformation parameters
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  integrand_return_t sector_`sectorID'_order_`cppOrder'_contour_deformation_polynomial#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  (#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const integration_variables,#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const complex_parameters,#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const deformation_parameters#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  )#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  {#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "integrand_return_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* c++ define the calls to "SecDecInternalRealPart"
  #redefine function "SecDecInternalRealPart"
  #Do labelID = 1, `largestLabel`function''
    Bracket SecDecInternalLabel`function'Call`labelID';
    .sort
    L realPart = expressionF[SecDecInternalLabel`function'Call`labelID'];
    .sort
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternal`function'Call`labelID' = SecDecInternalRealPart(%%E);" realPart(#@no_split_expression@#)
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
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
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalSecDecInternalDeformed`IV'Call = %%e" deformedIV(#@no_split_expression@#)
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"
      multiply replace_(SecDecInternalLabelSecDecInternalDeformed`IV', 0);
      .sort

    #EndIf
  #EndDo
  drop deformedIV;

* write `SecDecInternalContourDeformationPolynomial'
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e" expressionF(#@no_split_expression@#)
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "return tmp;#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`realParameters',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`complexParameters',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`occurringDeformationParameters',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalRealPart#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the function and the namespace in the c++ file
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* Delete "expressionF" since it is no longer needed.
  #clearoptimize
  drop expressionF;
  .sort


* The code above implements the contour deformation. Below, we implement an
* algorithm that finds reasonable `deformation_parameters`.


* We are writing the c++ file "optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp"
* and the corresponding header "optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp".
* The header can already completely be written here:
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#ifndef `name'_codegen_optimize_deformation_parameters_sector_`sectorID'_order_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#define `name'_codegen_optimize_deformation_parameters_sector_`sectorID'_order_`cppOrder'_hpp_included#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#include \"`name'.hpp\"#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#include \"functions.hpp\"#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#include <cmath>#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#include <limits>#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#include <vector>#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "{#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "  secdecutil::SectorContainerWithDeformation<real_t, complex_t>::MaximalDeformationFunction sector_`sectorID'_order_`cppOrder'_maximal_allowed_deformation_parameters;#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "};#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp> "#endif#@SecDecInternalNewline@#"

* Open the namespace and the function in the corresponding .cpp file
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#include \"optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.hpp\"#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "{#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "  void sector_`sectorID'_order_`cppOrder'_maximal_allowed_deformation_parameters#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "  (#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "      real_t * output_deformation_parameters,#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const integration_variables,#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "      complex_t const * const complex_parameters#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "  )#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "  {#@SecDecInternalNewline@#"

* We try to maximize the contour deformation by demandaing that the ratio
* of the second derivatives (the entries in the Hessian matrix) over the
* first derivatives is less than one. Each entry in the Hessian divided by
* a first derivative gives an upper bound for the corresponding deformation
* parameter.

* Construct the gradient and the Hessian
  L derivatives = SecDecInternalsDUMMYDerivativesAppendix;
  #Do idx1 = {`occurringIntegrationVariableIndices',}
    #If x`idx1' != x
      Id SecDecInternalsDUMMYDerivativesAppendix = SecDecInternalLabelGradient^(`idx1'+1) * d`SecDecInternalContourDeformationPolynomial'd`idx1'(`integrationVariables',`regulators') + SecDecInternalsDUMMYDerivativesAppendix;
      #Do idx2 = {`occurringIntegrationVariableIndices',}
        #If x`idx2' != x
          #If `idx1' <= `idx2'
            Id SecDecInternalsDUMMYDerivativesAppendix = SecDecInternalLabelHessianI^(`idx1'+1) * SecDecInternalLabelHessianJ^(`idx2'+1) * dd`SecDecInternalContourDeformationPolynomial'd`idx1'd`idx2'(`integrationVariables',`regulators') +
                                                         SecDecInternalsDUMMYDerivativesAppendix;
          #EndIf
        #EndIf
      #EndDo
    #EndIf
  #EndDo

* Set the following variables to zero:
*   - SecDecInternalsDUMMYDerivativesAppendix
*   - all regulators
*   - absent integration variables (if any)
  multiply replace_(
                       SecDecInternalsDUMMYDerivativesAppendix,0,
                       `nullifyRegulators'
                       #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
                         ,`nullifyAbsentIVs'
                       #EndIf
                   );
  .sort

* Explicitly insert the derivatives of `SecDecInternalContourDeformationPolynomial'
  argument SecDecInternalRealPart;
    #call insertOther
  endArgument;
  #call insertOther

* translate sympy's imaginary unit to FORM's imaginary unit
  multiply replace_(I,i_);
  .sort

* Define the integration variables, the real parameters, the complex parameters,
* and the deformation parameters as c preprocessor variables.
* (The c function takes them packed into an array).
* "Format rational": Need the indices as integers.
  Format rational;
  #call cppDefine(`occurringIntegrationVariables',integration_variables,optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`realParameters',real_parameters,optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`complexParameters',complex_parameters,optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)

* optimize
  Bracket SecDecInternalLabelGradient, SecDecInternalLabelHessianI, SecDecInternalLabelHessianJ;
  Format O`optimizationLevel';
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviation;
  #optimize derivatives

* Since FORM does not use "abbreviation[0]", we can use it as temporary variable.
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

* write the optimization symbols
  Format float 20;
  Format C;
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* define the entries of the gradient and the Hessian in c++
  Format rational;
  Format C;
  #$cppIdx1 = -1;
  #Do idx1 = {`occurringIntegrationVariableIndices',}
    #If x`idx1' != x
      #$cppIdx1 = $cppIdx1 + 1;
      Bracket SecDecInternalLabelGradient, SecDecInternalLabelHessianI, SecDecInternalLabelHessianJ;
      .sort
      L expr = derivatives[SecDecInternalLabelGradient^(`idx1'+1)];
      .sort
      #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
      #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_t abs_gradient_`$cppIdx1' = std::abs(tmp);#@SecDecInternalNewline@#"
      #$cppIdx2 = -1;
      #Do idx2 = {`occurringIntegrationVariableIndices',}
        #If x`idx2' != x
          #$cppIdx2 = $cppIdx2 + 1;
          #If `idx1' <= `idx2'
            Bracket SecDecInternalLabelGradient, SecDecInternalLabelHessianI, SecDecInternalLabelHessianJ;
            .sort
            L expr = derivatives[SecDecInternalLabelHessianI^(`idx1'+1)*SecDecInternalLabelHessianJ^(`idx2'+1)];
            .sort
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_t abs_Hessian_`$cppIdx1'_`$cppIdx2' = std::abs(tmp);#@SecDecInternalNewline@#"
          #EndIf
        #EndIf
      #EndDo
    #EndIf
  #EndDo

* c++-define the smallest nonzero number and a real temporary variable
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_t big = std::numeric_limits<real_t>::max();#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_t real_tmp;#@SecDecInternalNewline@#"

* find the maximal lambda that forces the second order to be smaller than
* the first order
* {

  #If `numOccurringIVOrder`shiftedOrderIndex'' != 0

    #$idx1 = -1;
    #Do IV1 = {`occurringIntegrationVariables'}
      #$idx1 = $idx1 + 1;

        #$idx2 = -1;
        #Do IV2 = {`occurringIntegrationVariables'}
          #$idx2 = $idx2 + 1;
          #$idx2plus1 = $idx2 + 1;

          #If `$idx1' <= `$idx2'
            #redefine HessianIDX "`$idx1'_`$idx2'"
          #Else
            #redefine HessianIDX "`$idx2'_`$idx1'"
          #EndIf

          #redefine nextDeformationParameterValue "abs_Hessian_`HessianIDX' == 0. ? big : abs_gradient_`$idx2' / abs_Hessian_`HessianIDX'"

          #If `$idx2' == 0
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "output_deformation_parameters[`$idx1'] = "
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "`nextDeformationParameterValue';#@SecDecInternalNewline@#"
          #Else
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_tmp = `nextDeformationParameterValue';#@SecDecInternalNewline@#"
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "output_deformation_parameters[`$idx1'] = "
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_tmp < output_deformation_parameters[`$idx1'] ?"
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_tmp : output_deformation_parameters[`$idx1'];#@SecDecInternalNewline@#"
          #EndIf

          #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

        #EndDo

    #EndDo

  #EndIf

* }

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables',optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`realParameters',optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`complexParameters',optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the function and the namespace in the c++ file
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* Delete the expression "derivatives" since it is no longer needed.
  #clearoptimize
  drop derivatives;
  .sort

* }
