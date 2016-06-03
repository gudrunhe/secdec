#-
Off statistics;

#include contour_deformation_sector`sectorID'.h

* We are writing a the c++ file "contour_deformation_sector_`sectorID'.cpp"
* and the corresponding header "contour_deformation_sector_`sectorID'.hpp".
* The header can already completely be written here:
#write <contour_deformation_sector_`sectorID'.hpp> "#ifndef __SecDec_include_guard_`name'_sector_`sectorID'_contourdef"
#write <contour_deformation_sector_`sectorID'.hpp> "#define __SecDec_include_guard_`name'_sector_`sectorID'_contourdef"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <`name'/config.hpp>"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <gsl/gsl_complex_math.h>"
#write <contour_deformation_sector_`sectorID'.hpp> "#include <gsl/gsl_linalg.h>"
#write <contour_deformation_sector_`sectorID'.hpp> "namespace `name'"
#write <contour_deformation_sector_`sectorID'.hpp> "{"
#write <contour_deformation_sector_`sectorID'.hpp> "  ContourDeformationFunction sector_`sectorID'__contourdef;"
#write <contour_deformation_sector_`sectorID'.hpp> "  IntegrandFunction sector_`sectorID'__F;"
#write <contour_deformation_sector_`sectorID'.hpp> "};"
#write <contour_deformation_sector_`sectorID'.hpp> "#endif"

* Open the namespace and the function in the corresponding .cpp file
#write <contour_deformation_sector_`sectorID'.cpp> "#include <`name'/integrands/contour_deformation_sector_`sectorID'.hpp>"
#write <contour_deformation_sector_`sectorID'.cpp> "namespace `name'"
#write <contour_deformation_sector_`sectorID'.cpp> "{"
#write <contour_deformation_sector_`sectorID'.cpp> "  integral_transformation_t sector_`sectorID'__contourdef"
#write <contour_deformation_sector_`sectorID'.cpp> "  ("
#write <contour_deformation_sector_`sectorID'.cpp> "      integration_variable_t const * const integration_variables,"
#write <contour_deformation_sector_`sectorID'.cpp> "      Mandelstam_t const * const Mandelstam,"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_mass_t const * const mass,"
#write <contour_deformation_sector_`sectorID'.cpp> "      deformation_parameter_t const * const deformation_parameters"
#write <contour_deformation_sector_`sectorID'.cpp> "  )"
#write <contour_deformation_sector_`sectorID'.cpp> "  {"

* The following procedure creates a local expression
* with the integral transformation (the "contour deformation")
* and its Jacobian matrix suitable for simultaneous optimization.
#call defineLocalExpression(contourdef)

* Implement commuting derivatives for F;
* i.e. replace ddFdjdi --> ddFdidj for i<=j.
* This is neccessary because pySecDec's algebra module does not
* imply commuting derivatives for general functions.
#$NumFPMinusOne = `NumFP' - 1;
.sort
#Do i = 0,`$NumFPMinusOne'
  #Do j = `i',`$NumFPMinusOne'
    Id ddFd`j'd`i'(?args) = ddFd`i'd`j'(?args);
  #EndDo
#EndDo

* Insert functions
Repeat;
  #call insert
EndRepeat;
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
#Do FP = {`FeynmanParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `FP' integration_variables[`$counter']"
  #$counter = $counter + 1;
#EndDo
#$counter = 0;
#Do DP = {`deformationParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `DP' deformation_parameters[`$counter']"
  #$counter = $counter + 1;
#EndDo
#$counter = 0;
#Do MS = {`Mandelstams'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `MS' Mandelstam[`$counter']"
  #$counter = $counter + 1;
#EndDo
#$counter = 0;
#Do m = {`masses'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#define `m' mass[`$counter']"
  #$counter = $counter + 1;
#EndDo

* Define "abbreviation[0]" as c++ preprocessor variable tmp.
* Since FORM does not use "abbreviation[0]", we can use it.
#write <contour_deformation_sector_`sectorID'.cpp> "#define tmp abbreviation[0]"

* write the transformation
* TODO: separate the Jacian matrix and the transformation
* TODO: calculate the Jacobian determinant numerically
Format float 20;
Format C;
Bracket labelTransform, labelJacobiI, labelJacobiJ;
.sort
#write <contour_deformation_sector_`sectorID'.cpp> "Feynman_t abbreviation[`optimmaxvar_' + 1];"
#write <contour_deformation_sector_`sectorID'.cpp> "%%O"
#write <contour_deformation_sector_`sectorID'.cpp> ""
#write <contour_deformation_sector_`sectorID'.cpp> "std::vector<Feynman_t> transformed_parameters(`NumFP');"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex *Jacobian = "
#write <contour_deformation_sector_`sectorID'.cpp> "    gsl_matrix_complex_alloc(`NumFP',`NumFP');"
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
    #write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex_set"
    #write <contour_deformation_sector_`sectorID'.cpp> "("
    #write <contour_deformation_sector_`sectorID'.cpp> "    Jacobian,"
    #write <contour_deformation_sector_`sectorID'.cpp> "    `i'-1, `j'-1,"
    #write <contour_deformation_sector_`sectorID'.cpp> "    {tmp.real(), tmp.imag()}"
    #write <contour_deformation_sector_`sectorID'.cpp> ");"
  #EndDo
#EndDo
* calculate the determinant numerically using the gsl
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_permutation * permutation = gsl_permutation_alloc (`NumFP');"
#write <contour_deformation_sector_`sectorID'.cpp> "int signum;"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_linalg_complex_LU_decomp (Jacobian, permutation, &signum);"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_complex det = gsl_linalg_complex_LU_det(Jacobian, signum);"
#write <contour_deformation_sector_`sectorID'.cpp> "tmp = {GSL_REAL(det), GSL_IMAG(det)};"

* free manually allocated memory
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_permutation_free(permutation);"
#write <contour_deformation_sector_`sectorID'.cpp> "gsl_matrix_complex_free(Jacobian);"

#write <contour_deformation_sector_`sectorID'.cpp> "return integral_transformation_t"
#write <contour_deformation_sector_`sectorID'.cpp> "    {"
#write <contour_deformation_sector_`sectorID'.cpp> "        std::move(transformed_parameters),"
#write <contour_deformation_sector_`sectorID'.cpp> "        std::move(tmp)"
#write <contour_deformation_sector_`sectorID'.cpp> "    };"

* Close the function in the c++ file
#write <contour_deformation_sector_`sectorID'.cpp> "  };"

* Optmimize the 2nd Symnanzik polynomial "F"
#clearoptimize
drop contourdef;
L expressionF = `F';
.sort
#optimize expressionF

* Write the function to optimize the contour deformation parameters
#write <contour_deformation_sector_`sectorID'.cpp> "  integrand_return_t sector_`sectorID'__F"
#write <contour_deformation_sector_`sectorID'.cpp> "  ("
#write <contour_deformation_sector_`sectorID'.cpp> "      Feynman_t const * const integration_variables,"
#write <contour_deformation_sector_`sectorID'.cpp> "      Mandelstam_t const * const Mandelstam,"
#write <contour_deformation_sector_`sectorID'.cpp> "      real_mass_t const * const mass"
#write <contour_deformation_sector_`sectorID'.cpp> "  )"
#write <contour_deformation_sector_`sectorID'.cpp> "  {"
#write <contour_deformation_sector_`sectorID'.cpp> "integrand_return_t abbreviation[`optimmaxvar_'];"
#write <contour_deformation_sector_`sectorID'.cpp> "%%O"
#write <contour_deformation_sector_`sectorID'.cpp> ""
#write <contour_deformation_sector_`sectorID'.cpp> "tmp = %%e" expressionF(tmp)
#write <contour_deformation_sector_`sectorID'.cpp> "return tmp;"

* undefine the c preprocessor macros
#Do FP = {`FeynmanParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `FP'"
#EndDo
#Do DP = {`deformationParameters'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `DP'"
#EndDo
#Do MS = {`Mandelstams'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `MS'"
#EndDo
#Do m = {`masses'}
  #write <contour_deformation_sector_`sectorID'.cpp> "#undef `m'"
#EndDo
#write <contour_deformation_sector_`sectorID'.cpp> "#undef tmp"

* Close the function and the namespace in the c++ file
#write <contour_deformation_sector_`sectorID'.cpp> "  };"
#write <contour_deformation_sector_`sectorID'.cpp> "};"

.end
