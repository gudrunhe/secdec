* The "lambda" parameters controlling the size of the deformation
#define deformationParameters "%(deformation_parameters)s"
Symbols `deformationParameters';

* The deformed integration variable functions (including appearing derivatives)
#define deformedIntegrationVariableDerivativeFunctions "%(deformed_integration_variable_derivative_functions)s"
CFunctions `deformedIntegrationVariableDerivativeFunctions';

* The Jacobian determinant of the contour deformation (including appearing derivatives)
#define contourdefJacobianFunctions "%(contourdef_Jacobian_derivative_functions)s"
CFunctions `contourdefJacobianFunctions';

* Define the calls to the contour deformation.
#Do function = {`deformedIntegrationVariableDerivativeFunctions',`contourdefJacobianFunctions'}
  AutoDeclare Symbols SecDecInternal`function'Call;
#EndDo

* Define the function that takes the real part
CFunction SecDecInternalRealPart;

* Define the call replacement symbols for the real part
AutoDeclare Symbols SecDecInternalSecDecInternalRealPartCall;

* Define the name of the polynomial for the contour deformation
* ("F" in loop integrals)
#define SecDecInternalContourDeformationPolynomial "%(contour_deformation_polynomial)s"

* Define the polynomials that should remain positive
* (e.g. "U" in loop integrals)
#define positivePolynomials "%(positive_polynomials)s"

* The transformation of the Feynman parameters
#procedure insertDeformedIntegrationVariables
  %(insert_deformed_integration_variables_procedure)s
#endProcedure

* Procedure that inserts the Jacobian determinant and
* its required derivatives. This procedure is written
* by python.
#procedure insertContourdefJacobianDerivatives
  %(insert_contourdef_Jacobian_derivatives_procedure)s
#endProcedure

* Procedure that removes vanishing derivatives of the deformed
* integration variables. This procedure is written by python.
#procedure removeVanishingDeformedIntegrationVariableDerivatives
  %(nullify_vanishing_deformed_integration_variable_calls_procedure)s
#endProcedure
