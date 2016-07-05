#-
Off statistics;

#include contour_deformation_sector`sectorID'.h

* We are writing a the c++ file "contour_deformation_sector_`sectorID'.cpp"
* and the corresponding header "contour_deformation_sector_`sectorID'.hpp".
* The header can already completely be written here:
#write <contour_deformation_sector_`sectorID'.hpp> "#ifndef `name'_codegen_contour_deformation_sector_`sectorID'_hpp_included#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#define `name'_codegen_contour_deformation_sector_`sectorID'_hpp_included#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <`name'/config.hpp>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <gsl/gsl_complex_math.h>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <gsl/gsl_linalg.h>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "namespace `name'#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "{#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "  ContourDeformationFunction contour_deformation;#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "  DeformableIntegrandFunction contour_deformation_polynomial;#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "};#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.hpp> "#endif#@SecDecInternalNewline@#"

* Open the namespace and the function in the corresponding .cpp file
#write <contour_deformation_sector_`sectorID'.cpp> "#include <`name'/integrands/contour_deformation_sector_`sectorID'.hpp>#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "namespace `name'#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "{#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  integral_transformation_t contour_deformation#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  (#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const integration_variables,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const deformation_parameters#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_t const * const real_parameters,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      complex_t const * const complex_parameters,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  )#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  {#@SecDecInternalNewline@#"

* The following procedure creates a local expression
* with the integral transformation (the "contour deformation")
* and its Jacobian matrix suitable for simultaneous optimization.
#call defineLocalExpression(contourdef)

* Implement commuting derivatives for the contour deformation polynomial F;
* i.e. replace ddFdjdi --> ddFdidj for i<=j.
* This is neccessary because pySecDec's algebra module does not
* imply commuting derivatives for general functions.
#$NumFPMinusOne = `NumFP' - 1;
.sort
#Do i = 0,`$NumFPMinusOne'
  #Do j = `i',`$NumFPMinusOne'
    Id dd`F'd`j'd`i'(?args) = dd`F'd`i'd`j'(?args);
  #EndDo
#EndDo

* define the argument for "replace_" that sets all regulators to zero
#redefine replaceArg ""
#$counter = 1;
#Do var = {`regulators',}
  #If `$counter' != 1
    #redefine replaceArg "`replaceArg' , "
  #EndIf
  #redefine replaceArg "`replaceArg' `var',0"
  #$counter = $counter + 1;
#EndDo

* Insert derivatives
#Do i = 0,`$NumFPMinusOne'
  #call generateReplacementd`F'd`i'(`replaceArg')
  .sort
  drop replacement;
  Id d`F'd`i'(?args) = replacement;
  .sort
  #Do j = `i',`$NumFPMinusOne'
    #call generateReplacementdd`F'd`i'd`j'(`replaceArg')
    .sort
    drop replacement;
    Id dd`F'd`i'd`j'(?args) = replacement;
    .sort
  #EndDo
#EndDo
.sort

* Simultaneously optimize the transformation and all
* entries of its Jacobi matrix.
Bracket labelTransform, labelJacobiI, labelJacobiJ;
Format O`optimizationLevel';
.sort
#optimize contourdef

* Define the integration, Mandelstam, and mass symbols as c preprocessor variables
* (The c function takes them packed into an array).
* "Format rational": Need the indices as integers.
Format rational;
#$counter = 0;
#Do IV = {`integrationVariables'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `IV' integration_variables[`$counter']#@SecDecInternalNewline@#"
  #$counter = $counter + 1;
#EndDo
#$counter = 0;
#Do DP = {`deformationParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `DP' deformation_parameters[`$counter']#@SecDecInternalNewline@#"
  #$counter = $counter + 1;
#EndDo
#$counter = 0;
#Do RP = {`realParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `RP' real_parameters[`$counter']#@SecDecInternalNewline@#"
  #$counter = $counter + 1;
#EndDo
#$counter = 0;
#Do CP = {`complexParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `CP' complex_parameters[`$counter']#@SecDecInternalNewline@#"
  #$counter = $counter + 1;
#EndDo

* Define "abbreviation[0]" as c++ preprocessor variable tmp.
* Since FORM does not use "abbreviation[0]", we can use it.
#write <contour_deformation_sector_`sectorID'.cpp> "#define tmp abbreviation[0]#@SecDecInternalNewline@#"

* write the transformation **** TODO: Only take the real part of `F' and its derivatives
Format float 20;
Format C;
Bracket labelTransform, labelJacobiI, labelJacobiJ;
.sort
#write <contour_deformation_sector_`sectorID'.cpp> "Feynman_t abbreviation[`optimmaxvar_' + 1];#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "%%O#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "std::vector<Feynman_t> transformed_parameters(`NumFP');#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex *Jacobian = #@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "    gsl_matrix_complex_alloc(`NumFP',`NumFP');#@SecDecInternalNewline@#"
#Do i = 1, `NumFP'
  L currentExpr = contourdef[labelTransform^`i'];
  Bracket labelTransform, labelJacobiI, labelJacobiJ;
  .sort
  #write <contour_deformation_sector_`sectorID'.cpp> "transformed_parameters[`i'-1] = %%e" currentExpr(transformed_parameters[`i'-1])
  #Do j = 1, `NumFP'
    L currentExpr = contourdef[labelJacobiI^`i' * labelJacobiJ^`j'];
    Bracket labelTransform, labelJacobiI, labelJacobiJ;
    .sort
    #write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e" currentExpr(tmp)
    #write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex_set#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "(#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "    Jacobian,#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "    `i'-1, `j'-1,#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> "    {tmp.real(), tmp.imag()}#@SecDecInternalNewline@#"
    #write <contour_deformation_sector_`sectorID'.cpp> ");#@SecDecInternalNewline@#"
  #EndDo
#EndDo
* calculate the determinant numerically using the gsl
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_permutation * permutation = gsl_permutation_alloc (`NumIV');#@SecDecInternalNewline@#"
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

* Optmimize the 2nd Symnanzik polynomial "F#@SecDecInternalNewline@#"
#clearoptimize
drop contourdef;
L expressionF = `F';
.sort
#optimize expressionF

* Write the function to optimize the contour deformation parameters
#write <contour_deformation_sector_`sectorID'.cpp> "  integrand_return_t contour_deformation_polynomial#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  (#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      Feynman_t const * const integration_variables,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      Mandelstam_t const * const Mandelstam,#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_mass_t const * const mass#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  )#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "  {#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "integrand_return_t abbreviation[`optimmaxvar_'];#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "%%O#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e" expressionF(tmp)
#write <contour_deformation_sector_`sectorID'.cpp> "return tmp;#@SecDecInternalNewline@#"

* undefine the c preprocessor macros
#Do FP = {`FeynmanParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `FP'#@SecDecInternalNewline@#"
#EndDo
#Do DP = {`deformationParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `DP'#@SecDecInternalNewline@#"
#EndDo
#Do MS = {`Mandelstams'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `MS'#@SecDecInternalNewline@#"
#EndDo
#Do m = {`masses'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `m'#@SecDecInternalNewline@#"
#EndDo
#write <contour_deformation_sector_`sectorID'.cpp> "#undef tmp#@SecDecInternalNewline@#"

* Close the function and the namespace in the c++ file
#write <contour_deformation_sector_`sectorID'.cpp> "  };#@SecDecInternalNewline@#"
#write <contour_deformation_sector_`sectorID'.cpp> "};#@SecDecInternalNewline@#"

.end
