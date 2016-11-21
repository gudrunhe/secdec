* The "lambda" parameters controlling the size of the deformation
#define deformationParameters "%(deformation_parameters)s"
Symbols `deformationParameters';

* The deformed integration variables
#define deformedIVNames "%(deformed_integration_variable_names)s"
CFunctions `deformedIVNames';

* The additional parameter controlling how fast the contour deformation goes to
* zero at the end points
Symbol SecDecInternalMu;

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

* We typically use commutativity of derivative to sort them.
* However, for contour deformation we need the unsorted second
* derivatives of the "SecDecInternalContourDeformationPolynomial",
* e.g. "ddFd1d0" in addition to "ddFd0d1".
#$numIVMinusOne = `numIV' - 1;
#Do i = 0,`$numIVMinusOne'
  #Do j = 0,`$numIVMinusOne'
    CFunctions dd`SecDecInternalContourDeformationPolynomial'd`i'd`j';
    Symbol SecDecInternaldd`SecDecInternalContourDeformationPolynomial'd`i'd`j'Call;
  #EndDo
#EndDo

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
