* The "lambda" parameters controlling the size of the deformation
#define deformationParameters "%(deformation_parameters)s"
Symbols `deformationParameters';

* The deformed integration variable functions (including appearing derivatives)
#define deformedIVFunctions "%(deformed_integration_variable_functions)s"
CFunctions `deformedIVFunctions';

* The additional parameter controlling how fast the contour deformation goes to
* zero at the end points
Symbol SecDecInternalMu;

* We multiply the error token by terms that we expect to vanish
Symbol SecDecInternalErrorToken;

* Define the function that takes the real part
CFunction SecDecInternalRealPart;

* Define the function call to the Jacobian determinant
CFunction SecDecInternalContourdefJacobian;

* Define the calls to the contour deformation.
#Do function = {`integrationVariables'}
  AutoDeclare Symbols SecDecInternalSecDecInternalDeformed`function'Call;
#EndDo
AutoDeclare Symbols SecDecInternalSecDecInternalContourdefJacobianCall;

* Define the function appearing in the contour deformation
CFunctions  SecDecIternalRealPart,
            SecDecInternalExpMinusMuOverX,
            SecDecInternalXExpMinusMuOverX,
           dSecDecInternalXExpMinusMuOverXd1;

* Define the call replacement symbols for the real part
AutoDeclare Symbols SecDecInternalSecDecInternalRealPartCall;

* Define the name of the polynomial for the contour deformation
* ("F" in loop integrals)
#define SecDecInternalContourDeformationPolynomial "%(contour_deformation_polynomial)s"

* The transformation of the Feynman parameters
#procedure insertDeformedIntegrationVariables
  %(insert_deformed_integration_variables_procedure)s
#endProcedure

* Procedure that inserts the Jacobian matrix
* of the contour deformation. This procedure
* is written by python.
#procedure insertContourdefJacobianMatrix
  %(insert_contourdef_Jacobian_procedure)s
#endProcedure

* Procedure that inserts the derivatives of the
* Jacobian matrix of the contour deformation. This
* procedure is written by python.
#procedure insertContourdefJacobianDerivatives
  %(insert_contourdef_Jacobian_derivatives_procedure)s
#endProcedure
