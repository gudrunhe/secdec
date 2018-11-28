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
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#ifdef SECDEC_WITH_CUDA#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) (complex_t{x}).real()#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#else#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) std::real(x)#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#endif#@SecDecInternalNewline@#"

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
  #call cppDefine(`occurringIntegrationVariables',integration_variables,optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`realParameters',real_parameters,optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`complexParameters',complex_parameters,optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#ifdef SECDEC_WITH_CUDA#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) (complex_t{x}).real()#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalAbs(x) thrust::abs(complex_t{x})#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#else#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalRealPart(x) std::real(x)#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#define SecDecInternalAbs(x) std::abs(x)#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#endif#@SecDecInternalNewline@#"

* optimize
  Bracket SecDecInternalLabelDeformation;
  Format O`optimizationLevel';
  .sort
  ExtraSymbols,array,SecDecInternalAbbreviation;
  #optimize deformations
  intohide deformations;
  Bracket SecDecInternalLabelDeformation;
  .sort

* Since FORM does not use "abbreviation[0]", we can use it as temporary variable.
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

* write the optimization symbols
  Format float 20;
  Format C;
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* set the output lambdas to ``1 / Abs(    Re(  x_k * (1-x_k) * dF_dx_k  )    )``
  #$cppidx = -1;
  #Do idx = {`occurringIntegrationVariableIndices',}
    #If x`idx' != x
      #$cppidx = $cppidx + 1;
      L expr = deformations[SecDecInternalLabelDeformation^(`idx'+1)];
      .sort
      Format rational;
      Format C;
      #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "output_deformation_parameters[`$cppidx'] ="
      Format float 20;
      Format C;
      #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "1./SecDecInternalAbs(SecDecInternalRealPart(%%E));#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
    #EndIf
  #EndDo

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables',optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`realParameters',optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`complexParameters',optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp)
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalRealPart#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalAbs#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the function and the namespace in the c++ file
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* Delete the expressions "expr" and "deformations" since they are no longer needed.
  #clearoptimize
  unhide deformations;
  drop expr;
  drop deformations;
  .sort

* }
