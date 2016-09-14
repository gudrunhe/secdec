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
CFunctions  SecDecInternalExpMinusMuOverX,
            SecDecInternalXExpMinusMuOverX,
           dSecDecInternalXExpMinusMuOverXd1;

* Define the name of the polynomial for the contour deformation
* ("F" in loop integrals)
#define F "%(contour_deformation_polynomial)s"

* We typically use commutativity of derivative to sort them.
* However, for contour deformation we need the unsorted second
* derivatives of `F', e.g. "ddFd1d0" in addition to "ddFd0d1".
#$numIVMinusOne = `numIV' - 1;
#Do i = 0,`$numIVMinusOne'
  #Do j = 0,`$numIVMinusOne'
    CFunctions dd`F'd`i'd`j', SecDecInternaldd`F'd`i'd`j';
    Symbol SecDecInternaldd`F'd`i'd`j'Call;
  #EndDo
#EndDo

* The transformation of the Feynman parameters
#procedure insertDeformedIntegrationVariables
  %(insert_deformed_integration_variables_procedure)s
#endProcedure

* The transformation of the Feynman parameters
* and its Jacobian matrix suitable for simultaneous
* code optimization. This expression is written by python.
#procedure insertContourdefExpression
  %(contourdef_expression_definition_procedure)s
#endProcedure
