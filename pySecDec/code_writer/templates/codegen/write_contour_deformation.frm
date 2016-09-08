* This file is included in "write_integrand.frm".
* It should be seen as a procedure that writes
* the code for the contour deformation.
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
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_complex_math.h>#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#include <gsl/gsl_linalg.h>#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "namespace `name'#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "{#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "  secdecutil::SectorContainerWithDeformation<real_t, complex_t>::ContourDeformationFunction sector_`sectorID'_order_`cppOrder'_contour_deformation;#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "  secdecutil::SectorContainerWithDeformation<real_t, complex_t>::DeformableIntegrandFunction sector_`sectorID'_order_`cppOrder'_contour_deformation_polynomial;#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "};#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.hpp> "#endif#@SecDecInternalNewline@#"

* Open the namespace and the function in the corresponding .cpp file
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

* The following procedure creates a local expression
* with the integral transformation (the "contour deformation")
* and its Jacobian matrix suitable for simultaneous optimization.
  Local contourdef = SecDecInternalsDUMMYContourdefExpression + SecDecInternalsDUMMYContourdefAppendix;
  #call insertContourdefExpression
  .sort

* Implement commuting derivatives for the contour deformation polynomial F;
* i.e. replace ddFdjdi --> ddFdidj for i<=j.
* This is neccessary because pySecDec's algebra module does not
* imply commuting derivatives for general functions.
  #Do i = 0,`$numIVMinusOne'
    #Do j = `i',`$numIVMinusOne'
      Id dd`F'd`j'd`i'(?args) = dd`F'd`i'd`j'(?args);
    #EndDo
  #EndDo
  .sort

* We must take only the real part of `F' and its derivatives.
* Collect all these and replace calls by symbols.
  Id `F'(?args) = SecDecInternal`F'Call;
  Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                          `F'(`integrationVariables',`regulators')*SecDecInternalLabel`F';
  #Do idx1 = {`occurringIntegrationVariableIndices',}
    #If x`idx1' != x
      Id d`F'd`idx1'(?args) = SecDecInternald`F'd`idx1'Call;
      Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                            d`F'd`idx1'(`integrationVariables',`regulators')*SecDecInternalLabeld`F'd`idx1';
      #Do idx2 = {`occurringIntegrationVariableIndices',}
        #If x`idx2' != x
          #If `idx1' <= `idx2'
            Id dd`F'd`idx1'd`idx2'(?args) = SecDecInternaldd`F'd`idx1'd`idx2'Call;
            Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                                    dd`F'd`idx1'd`idx2'(`integrationVariables',`regulators')*SecDecInternalLabeldd`F'd`idx1'd`idx2';
          #EndIf
        #EndIf
      #EndDo
    #EndIf
  #EndDo
  Id SecDecInternalsDUMMYContourdefAppendix = 0;

* Remove derivatives of `F' by absent integration variables
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

* set all regulators and "SecDecInternalsDUMMYToOptimize" to zero
  multiply replace_(`nullifyRegulators' , SecDecInternalsDUMMYToOptimize,0);
  .sort

  Format rational;
  Format Normal;

* define the argument for "replace_" that sets all absent integration variables to zero
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

* Remove absent integration variables
  #If `numIV' != `numOccurringIVOrder`shiftedOrderIndex''
    multiply replace_(`nullifyAbsentIVs');
    .sort
  #EndIf

* Explicitly insert `F' and its derivatives
  #call insert

* translate sympy's imaginary unit to FORM's imaginary unit
  multiply replace_(I,i_);
  .sort

* Define the integration variables, the real parameters, the complex parameters,
* and the deformation parameters as c preprocessor variables.
* (The c function takes them packed into an array).
* "Format rational": Need the indices as integers.
  Format rational;
  #call cppDefine(`occurringIntegrationVariables',integration_variables,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`realParameters',real_parameters,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`complexParameters',complex_parameters,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppDefine(`occurringDeformationParameters',deformation_parameters,contour_deformation_sector_`sectorID'_`cppOrder'.cpp)

* optimize
  AB `integrationVariables', `realParameters', `complexParameters';
  Format O`optimizationLevel';
  .sort
  #optimize contourdef

* Since FORM does not use "abbreviation[0]", we can use it as temporary variable.
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

* write the optimization symbols
  Format float 20;
  Format C;
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "complex_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "%%O#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#@SecDecInternalNewline@#"

* define the symbols "SecDecInternal...Call" as c++ variables
* {

  Format rational;
  Format C;
  Bracket SecDecInternalLabel`F';
  .sort
  L expr = contourdef[SecDecInternalLabel`F'];
  .sort
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternal`F'Call = tmp.real();#@SecDecInternalNewline@#"

  #Do idx1 = {`occurringIntegrationVariableIndices',}
    #If x`idx1' != x
      Bracket SecDecInternalLabeld`F'd`idx1';
      .sort
      L expr = contourdef[SecDecInternalLabeld`F'd`idx1'];
      .sort
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
      #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternald`F'd`idx1'Call = tmp.real();#@SecDecInternalNewline@#"

      #Do idx2 = {`occurringIntegrationVariableIndices',}
        #If x`idx2' != x
          #If `idx1' <= `idx2'
            Bracket SecDecInternalLabeldd`F'd`idx1'd`idx2';
            .sort
            L expr = contourdef[SecDecInternalLabeldd`F'd`idx1'd`idx2'];
            .sort
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(#@no_split_expression@#)
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "real_t SecDecInternaldd`F'd`idx1'd`idx2'Call = tmp.real();#@SecDecInternalNewline@#"
          #EndIf
        #EndIf
      #endDo
    #EndIf
  #endDo

* }

* write the transformation
  Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
  .sort
  Format float 20;
  Format C;
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "std::vector<complex_t> transformed_parameters(`numOccurringIVOrder`shiftedOrderIndex'');#@SecDecInternalNewline@#"

  #If `numOccurringIVOrder`shiftedOrderIndex'' != 0

    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex *Jacobian = #@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    gsl_matrix_complex_alloc(`numOccurringIVOrder`shiftedOrderIndex'',`numOccurringIVOrder`shiftedOrderIndex'');#@SecDecInternalNewline@#"
    #$i = -1;
    #Do idx1 = {`occurringIntegrationVariableIndices',}
      #If x`idx1' != x
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
        #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "transformed_parameters[`$i'] = %%e" expr(#@no_split_expression@#)
        #$j = -1;
        #Do idx2 = {`occurringIntegrationVariableIndices',}
          #If x`idx2' != x
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
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e" expr(#@no_split_expression@#)
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex_set#@SecDecInternalNewline@#"
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "(#@SecDecInternalNewline@#"
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    Jacobian,#@SecDecInternalNewline@#"
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    `$i', `$j',#@SecDecInternalNewline@#"
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    {tmp.real(), tmp.imag()}#@SecDecInternalNewline@#"
            #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> ");#@SecDecInternalNewline@#"
          #EndIf
        #EndDo
      #EndIf
    #EndDo

*   calculate the determinant numerically using the gsl
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_permutation * permutation = gsl_permutation_alloc (`numOccurringIVOrder`shiftedOrderIndex'');#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "int signum;#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_linalg_complex_LU_decomp (Jacobian, permutation, &signum);#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_complex det = gsl_linalg_complex_LU_det(Jacobian, signum);#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = {GSL_REAL(det), GSL_IMAG(det)};#@SecDecInternalNewline@#"

*   free manually allocated memory
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_permutation_free(permutation);#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "gsl_matrix_complex_free(Jacobian);#@SecDecInternalNewline@#"

  #Else

    #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = 1;#@SecDecInternalNewline@#"

  #Endif

  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "return secdecutil::integral_transformation_t<complex_t>#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    {#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "        std::move(transformed_parameters),#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "        std::move(tmp)#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "    };#@SecDecInternalNewline@#"

* Close the function in the c++ file
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"

* Optmimize the 2nd Symnanzik polynomial "F"
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

* Write the function to optimize the contour deformation parameters
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
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "tmp = %%e" expressionF(#@no_split_expression@#)
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "return tmp;#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
  #call cppUndefine(`occurringIntegrationVariables',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`realParameters',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`complexParameters',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #call cppUndefine(`occurringDeformationParameters',contour_deformation_sector_`sectorID'_`cppOrder'.cpp)
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#undef SecDecInternalMu#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the function and the namespace in the c++ file
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "  };#@SecDecInternalNewline@#"
  #write <contour_deformation_sector_`sectorID'_`cppOrder'.cpp> "};#@SecDecInternalNewline@#"

* Delete the expression "contourdef" since it is no longer needed.
  drop contourdef;
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
      Id SecDecInternalsDUMMYDerivativesAppendix = SecDecInternalLabelGradient^(`idx1'+1) * d`F'd`idx1'(`integrationVariables',`regulators') + SecDecInternalsDUMMYDerivativesAppendix;
      #Do idx2 = {`occurringIntegrationVariableIndices',}
        #If x`idx2' != x
          #If `idx1' <= `idx2'
            Id SecDecInternalsDUMMYDerivativesAppendix = SecDecInternalLabelHessianI^(`idx1'+1) * SecDecInternalLabelHessianJ^(`idx2'+1) * dd`F'd`idx1'd`idx2'(`integrationVariables',`regulators') +
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

* Explicitly insert the derivatives of `F'
  #call insert

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
  #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_t tiny = std::numeric_limits<real_t>::min();#@SecDecInternalNewline@#"
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

          #If `$idx2' == 0
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "output_deformation_parameters[`$idx1'] = "
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "abs_gradient_`$idx2' / (abs_Hessian_`HessianIDX' + tiny);#@SecDecInternalNewline@#"
          #Else
            #write <optimize_deformation_parameters_sector_`sectorID'_`cppOrder'.cpp> "real_tmp = abs_gradient_`$idx2' / (abs_Hessian_`HessianIDX' + tiny);#@SecDecInternalNewline@#"
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
