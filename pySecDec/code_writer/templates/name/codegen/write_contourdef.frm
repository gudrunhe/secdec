#-
Off statistics;

* Define two general procedures that write c++ code to define and undefine
* c++ preprocessor varibables accessing a c++ array.
#procedure cppDefine(?FORMNames,cppArrayName)
  #$counter = 0;
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <contour_deformation_sector_`sectorID'.cpp> "#define `varname' `cppArrayName'[`$counter']#@SecDecInternalNewline@#"
      #$counter = $counter + 1;
    #EndIf
  #EndDo
#endProcedure
#procedure cppUndefine(?FORMNames)
  #Do varname = {`?FORMNames'}
    #If x`varname' != x
      #write <contour_deformation_sector_`sectorID'.cpp> "#undef `varname'#@SecDecInternalNewline@#"
    #EndIf
  #EndDo
#endProcedure

#include contour_deformation_sector`sectorID'.h

* We are writing the c++ file "contour_deformation_sector_`sectorID'.cpp"
* and the corresponding header "contour_deformation_sector_`sectorID'.hpp".
* The header can already completely be written here:
#write <contour_deformation_sector_`sectorID'.hpp> "#ifndef `name'_codegen_contour_deformation_sector_`sectorID'_hpp_included#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#define `name'_codegen_contour_deformation_sector_`sectorID'_hpp_included#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <`name'/config.hpp>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <gsl/gsl_complex_math.h>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <gsl/gsl_linalg.h>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "namespace `name'#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "{#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "  ContourDeformationFunction sector_`sectorID'_contour_deformation;#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "  DeformableIntegrandFunction sector_`sectorID'_contour_deformation_polynomial;#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "};#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#endif#@SecDecInternalNewline@#"

* Open the namespace and the function in the corresponding .cpp file
#write <contour_deformation_sector_`sectorID'.cpp> "#include <`name'/integrands/contour_deformation_sector_`sectorID'.hpp>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "namespace `name'#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "{#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  integral_transformation_t sector_`sectorID'_contour_deformation#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  (#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const integration_variables,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      complex_t const * const complex_parameters,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const deformation_parameters#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  )#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  {#@SecDecInternalNewline@#"

* The following procedure creates a local expression
* with the integral transformation (the "contour deformation")
* and its Jacobian matrix suitable for simultaneous optimization.
Local contourdef = SecDecInternalsDUMMYContourdefExpression + SecDecInternalsDUMMYContourdefAppendix;
#call insertContourdefExpression

* Implement commuting derivatives for the contour deformation polynomial F;
* i.e. replace ddFdjdi --> ddFdidj for i<=j.
* This is neccessary because pySecDec's algebra module does not
* imply commuting derivatives for general functions.
.sort
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
#Do i = 0,`$numIVMinusOne'
  Id d`F'd`i'(?args) = SecDecInternald`F'd`i'Call;
  Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                              d`F'd`i'(`integrationVariables',`regulators')*SecDecInternalLabeld`F'd`i';
  #Do j = 0,`i'
    Id dd`F'd`j'd`i'(?args) = SecDecInternaldd`F'd`j'd`i'Call;
    Id SecDecInternalsDUMMYContourdefAppendix = SecDecInternalsDUMMYContourdefAppendix +
                                                dd`F'd`j'd`i'(`integrationVariables',`regulators')*SecDecInternalLabeldd`F'd`j'd`i';
  #EndDo
#EndDo

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
#call cppDefine(`integrationVariables',integration_variables)
#call cppDefine(`realParameters',real_parameters)
#call cppDefine(`complexParameters',complex_parameters)
#call cppDefine(`deformationParameters',deformation_parameters)

* optimize
AB `integrationVariables', `realParameters', `complexParameters';
Format O`optimizationLevel';
.sort
#optimize contourdef

* Since FORM does not use "abbreviation[0]", we can use it as temporary variable.
#write <contour_deformation_sector_`sectorID'.cpp> "#define tmp SecDecInternalAbbreviation[0]#@SecDecInternalNewline@#"

* write the optimization symbols
Format float 20;
Format C;
#write <contour_deformation_sector_`sectorID'.cpp> "complex_t SecDecInternalAbbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "%%O#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"

* define the symbols "SecDecInternal...Call" as c++ variables
* {

Format rational;
Format C;
Bracket SecDecInternalLabel`F';
.sort
L expr = contourdef[SecDecInternalLabel`F'];
.sort
#write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(tmp)
#write <contour_deformation_sector_`sectorID'.cpp> "real_t SecDecInternal`F'Call = tmp.real();#@SecDecInternalNewline@#"

#Do i = 0,`$numIVMinusOne'
  Bracket SecDecInternalLabeld`F'd`i';
  .sort
  L expr = contourdef[SecDecInternalLabeld`F'd`i'];
  .sort
  #write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(tmp)
  #write <contour_deformation_sector_`sectorID'.cpp> "real_t SecDecInternald`F'd`i'Call = tmp.real();#@SecDecInternalNewline@#"

  #Do j = `i',`$numIVMinusOne'
    Bracket SecDecInternalLabeldd`F'd`i'd`j';
    .sort
    L expr = contourdef[SecDecInternalLabeldd`F'd`i'd`j'];
    .sort
    #write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e#@SecDecInternalNewline@#" expr(tmp)
    #write <contour_deformation_sector_`sectorID'.cpp> "real_t SecDecInternaldd`F'd`i'd`j'Call = tmp.real();#@SecDecInternalNewline@#"
  #endDo
#endDo

* }

* write the transformation
Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
.sort
Format float 20;
Format C;
#write <contour_deformation_sector_`sectorID'.cpp> "std::vector<complex_t> transformed_parameters(`numIV');#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex *Jacobian = #@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "    gsl_matrix_complex_alloc(`numIV',`numIV');#@SecDecInternalNewline@#"
#Do i = 1, `numIV'
  Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
  .sort
  L expr = contourdef[SecDecInternalLabelTransformation^`i'];
  .sort
  #write <contour_deformation_sector_`sectorID'.cpp> "transformed_parameters[`i'-1] = %%e" expr(transformed_parameters[`i'-1])
  #Do j = 1, `numIV'
    Bracket SecDecInternalLabelTransformation, SecDecInternalLabelJacobianMatrixI, SecDecInternalLabelJacobianMatrixJ;
    .sort
    L expr = contourdef[SecDecInternalLabelJacobianMatrixI^`i' * SecDecInternalLabelJacobianMatrixJ^`j'];
    .sort
    #write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e" expr(tmp)
    #write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex_set#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "(#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "    Jacobian,#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "    `i'-1, `j'-1,#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "    {tmp.real(), tmp.imag()}#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> ");#@SecDecInternalNewline@#"
  #EndDo
#EndDo
* calculate the determinant numerically using the gsl
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_permutation * permutation = gsl_permutation_alloc (`numIV');#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "int signum;#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_linalg_complex_LU_decomp (Jacobian, permutation, &signum);#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_complex det = gsl_linalg_complex_LU_det(Jacobian, signum);#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "tmp = {GSL_REAL(det), GSL_IMAG(det)};#@SecDecInternalNewline@#"

* free manually allocated memory
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_permutation_free(permutation);#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex_free(Jacobian);#@SecDecInternalNewline@#"

#write <contour_deformation_sector_`sectorID'.cpp> "return integral_transformation_t#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "    {#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "        std::move(transformed_parameters),#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "        std::move(tmp)#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "    };#@SecDecInternalNewline@#"

* Close the function in the c++ file
#write <contour_deformation_sector_`sectorID'.cpp> "  };#@SecDecInternalNewline@#"


* Optmimize the 2nd Symnanzik polynomial "F"
#clearoptimize
drop contourdef;
L expressionF = `F'(`integrationVariables',`regulators')*replace_(`nullifyRegulators');
.sort
#call insert
.sort
#optimize expressionF

* Write the function to optimize the contour deformation parameters
#write <contour_deformation_sector_`sectorID'.cpp> "  integrand_return_t sector_`sectorID'_contour_deformation_polynomial#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  (#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      complex_t const * const integration_variables,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      complex_t const * const complex_parameters#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  )#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  {#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "integrand_return_t SecDecInternalAbbreviation[`optimmaxvar_'];#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "%%O#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e" expressionF(tmp)
#write <contour_deformation_sector_`sectorID'.cpp> "return tmp;#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
#call cppUndefine(`integrationVariables')
#call cppUndefine(`realParameters')
#call cppUndefine(`complexParameters')
#call cppUndefine(`deformationParameters')
#write <contour_deformation_sector_`sectorID'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the function and the namespace in the c++ file
#write <contour_deformation_sector_`sectorID'.cpp> "  };#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "};#@SecDecInternalNewline@#"

.end
