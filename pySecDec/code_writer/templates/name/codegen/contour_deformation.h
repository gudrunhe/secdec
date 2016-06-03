#include sector`sectorID'.h

* We need labels for the code optimization
Symbols labelTransform, labelJacobiI, labelJacobiJ;

* The "lambda" parameters controlling the size of the deformation
#define deformationParameters "%(deformation_parameters)s"
Symbols `deformationParameters';

* The transformation of the Feynman parameters
* (z_k({x_k}) = x_k - i * lambda_k * (1 - x_k) * Re(dF_dx_k))
* and its Jacobian matrix suitable for simultaneous
* code optimization. This expression is written by python.
#procedure defineLocalExpression(contourdef)
  Local `contourdef' = %(contourdef_expression)s;
#endProcedure
