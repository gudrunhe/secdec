from pySecDec.code_writer.sum_package import Coefficient

### Coefficients definitions ###

parameters = ['s12','s13','s23','mtsq','mhsq']


#################################################################################################################################################################
####################################################################    FINAL EXPRESSIONS  ######################################################################
#################################################################################################################################################################

# P1  =  -(  (16*eps*(s12^2 - s13*s23)  )/  ( (1 - 2*eps)*(2 - 2*eps)*s12*(s12 + s13)*s23*(s12 + s23)  )    )  

P1num = ['-(16*eps*(s12^2 - s13*s23))']

P1den = ['((1 - 2*eps)*(2 - 2*eps)*s12*(s12 + s13)*s23*(s12 + s23))']

# P2  =  (8*eps*(4*mtsq - s12)) /  (  (1 - 2*eps)*(2 - 2*eps)*s12*s23  )  

P2num = ['(8*eps*(4*mtsq - s12))']

P2den = ['(  (1 - 2*eps)*(2 - 2*eps)*s12*s23  )']

# P3  =  -( (4*(4*mtsq - s13)*s13*(-2*eps*s12^2 + 2*(4 - 2*eps)*s12*s23 + (4 - 2*eps)*s23^2) )  / (  (1 - 2*eps)*(2 - 2*eps)*s12^2*s23*(s12 + s23)^2)  ) 

P3num = ['-(4*(4*mtsq - s13)*s13*(-2*eps*s12^2 + 2*(4 - 2*eps)*s12*s23 + (4 - 2*eps)*s23^2) ) ']

P3den = ['(  (1 - 2*eps)*(2 - 2*eps)*s12^2*s23*(s12 + s23)^2) ']

# P4  =   -((4*(-2*eps*s12^2 + 2*(4 - 2*eps)*s12*s13 + (4 - 2*eps)*s13^2)*(4*mtsq - s23))/  ((1 - 2*eps)*(2 - 2*eps)*s12^2*(s12 + s13)^2))  

P4num = ['-(4*(-2*eps*s12^2 + 2*(4 - 2*eps)*s12*s13 + (4 - 2*eps)*s13^2)*(4*mtsq - s23))']

P4den = ['((1 - 2*eps)*(2 - 2*eps)*s12^2*(s12 + s13)^2)']

# P6  = (2*eps*(s13^2 + s23^2))/ ((1 - 2*eps)*s13*s23^2)  

P6num = ['(2*eps*(s13^2 + s23^2))']

P6den = ['((1 - 2*eps)*s13*s23^2)']

# P7  =  -((s13^2*(-2*eps*s12^2 + (4 - 2*eps)*s23^2))/((1 - 2*eps)*s12^3*s23^2))   

P7num = ['-(s13^2*(-2*eps*s12^2 + (4 - 2*eps)*s23^2))']

P7den = ['((1 - 2*eps)*s12^3*s23^2)']

#  P8  =  -(((-2*eps*s12^2 + (4 - 2*eps)*s13^2)*s23)/((1 - 2*eps)*s12^3*s13)) 

P8num = ['-((-2*eps*s12^2 + (4 - 2*eps)*s13^2)*s23)']

P8den = ['((1 - 2*eps)*s12^3*s13)']

# P12 = (-4*mtsq*s13 + (1 - 2*eps)*s12*s13 - 2*eps*s12*s23)/((1 - 2*eps)*s12*s13)  

P12num = ['(-4*mtsq*s13 + (1 - 2*eps)*s12*s13 - 2*eps*s12*s23)']

P12den = ['((1 - 2*eps)*s12*s13)']

# P13 =  (s13*(-2*eps*s12*s13 - 4*mtsq*s23 + (1 - 2*eps)*s12*s23))/ ((1 - 2*eps)*s12*s23^2)

P13num = ['(s13*(-2*eps*s12*s13 - 4*mtsq*s23 + (1 - 2*eps)*s12*s23))']

P13den = ['((1 - 2*eps)*s12*s23^2)']

# P14 =  (s13*(12*mtsq*s12 - (1 - 2*eps)*s12^2 + (4 - 2*eps)*s13*s23))/  ((1 - 2*eps)*s12^3)

P14num = ['(s13*(12*mtsq*s12 - (1 - 2*eps)*s12^2 + (4 - 2*eps)*s13*s23))']

P14den = ['((1 - 2*eps)*s12^3)']

# P5 = (4*(-4*mtsq + s12 + s13 + s23)*(-2*s13*s23*(4*s12^3 + 5*s12^2*(s13 + s23) + s13*s23*(s13 + s23) + 2*s12*(s13 + s23)^2) + eps*(s12 + s13)*(s12 + s23)*(3*s12^3 + 2*s12^2*(s13 + s23) + s13*s23*(s13 + s23) + s12*(s13^2 + s13*s23 + s23^2))))  /   (  (-1 + eps)*(-1 + 2*eps)*s12^2*(s12 + s13)^2*s23*(s12 + s23)^2)

P5num = ['(4*(-4*mtsq + s12 + s13 + s23)*(-2*s13*s23*(4*s12^3 + 5*s12^2*(s13 + s23) + s13*s23*(s13 + s23) + 2*s12*(s13 + s23)^2) + eps*(s12 + s13)*(s12 + s23)*(3*s12^3 + 2*s12^2*(s13 + s23) + s13*s23*(s13 + s23) + s12*(s13^2 + s13*s23 + s23^2)))) ']

P5den = [' ((-1 + eps)*(-1 + 2*eps)*s12^2*(s12 + s13)^2*s23*(s12 + s23)^2)']

# P9 = (2*(s13 + s23)*(8*eps*mtsq*s13*s23 + (-4*mtsq + s12)*s13*s23 + eps^2*s12*(s13 + s23)^2 - eps*s12*(s13^2 + 3*s13*s23 + s23^2)))/  ( (-1 + eps)*(-1 + 2*eps)*s12^2*s13*s23^2)

P9num = ['(2*(s13 + s23)*(8*eps*mtsq*s13*s23 + (-4*mtsq + s12)*s13*s23 + eps^2*s12*(s13 + s23)^2 - eps*s12*(s13^2 + 3*s13*s23 + s23^2)))']

P9den = ['( (-1 + eps)*(-1 + 2*eps)*s12^2*s13*s23^2)']

# P10 = (2*(eps^2*(s12 + s23)*(s12^2*(s13 - 2*s23)*s23 + s12*s13*s23^2 + s13*s23^3 + s12^3*(s13 + 2*s23)) - eps*(2*s12^3*(-4*mtsq + s13)*s23 + s12^2*(16*mtsq + 4*s13 - 3*s23)*s23^2 + 2*s12*(4*mtsq +  3*s13)*s23^3 + 3*s13*s23^4 + s12^4*(s13 + 3*s23)) + s23*(4*mtsq*s12*(-s12^2 + 2*s12*s23 + s23^2) + (s12 + s23)*(s12^3 - s12^2*s23 + 2*s12*s13*s23 + 2*s13*s23^2))))/  ((-1 + eps)*(-1 + 2*eps)*s12^3*s23^2*(s12 + s23))

P10num = ['(2*(eps^2*(s12 + s23)*(s12^2*(s13 - 2*s23)*s23 + s12*s13*s23^2 + s13*s23^3 + s12^3*(s13 + 2*s23)) - eps*(2*s12^3*(-4*mtsq + s13)*s23 + s12^2*(16*mtsq + 4*s13 - 3*s23)*s23^2 + 2*s12*(4*mtsq +  3*s13)*s23^3 + 3*s13*s23^4 + s12^4*(s13 + 3*s23)) + s23*(4*mtsq*s12*(-s12^2 + 2*s12*s23 + s23^2) + (s12 + s23)*(s12^3 - s12^2*s23 + 2*s12*s13*s23 + 2*s13*s23^2))))']

P10den = ['((-1 + eps)*(-1 + 2*eps)*s12^3*s23^2*(s12 + s23))']

# P11 = (2*(-4*mtsq*s12^3*s13 + 8*eps*mtsq*s12^3*s13 + s12^4*s13 - 3*eps*s12^4*s13 + 2*eps^2*s12^4*s13 + 8*mtsq*s12^2*s13^2 - 16*eps*mtsq*s12^2*s13^2 + 4*mtsq*s12*s13^3 - 8*eps*mtsq*s12*s13^3 - s12^2*s13^3 +  3*eps*s12^2*s13^3 - 2*eps^2*s12^2*s13^3 - eps*s12^4*s23 +  eps^2*s12^4*s23 - 2*eps*s12^3*s13*s23 + 2*eps^2*s12^3*s13*s23 +   2*s12^2*s13^2*s23 - 4*eps*s12^2*s13^2*s23 +  2*eps^2*s12^2*s13^2*s23 +   4*s12*s13^3*s23 - 6*eps*s12*s13^3*s23 + 2*eps^2*s12*s13^3*s23 +   2*s13^4*s23 - 3*eps*s13^4*s23 + eps^2*s13^4*s23))/ ((-1 + eps)*(-1 + 2*eps)*s12^3*s13*(s12 + s13)*s23)

P11num = ['(2*(-4*mtsq*s12^3*s13 + 8*eps*mtsq*s12^3*s13 + s12^4*s13 - 3*eps*s12^4*s13 + 2*eps^2*s12^4*s13 + 8*mtsq*s12^2*s13^2 - 16*eps*mtsq*s12^2*s13^2 + 4*mtsq*s12*s13^3 - 8*eps*mtsq*s12*s13^3 - s12^2*s13^3 +  3*eps*s12^2*s13^3 - 2*eps^2*s12^2*s13^3 - eps*s12^4*s23 +  eps^2*s12^4*s23 - 2*eps*s12^3*s13*s23 + 2*eps^2*s12^3*s13*s23 +   2*s12^2*s13^2*s23 - 4*eps*s12^2*s13^2*s23 +  2*eps^2*s12^2*s13^2*s23 +   4*s12*s13^3*s23 - 6*eps*s12*s13^3*s23 + 2*eps^2*s12*s13^3*s23 +   2*s13^4*s23 - 3*eps*s13^4*s23 + eps^2*s13^4*s23))']

P11den = ['((-1 + eps)*(-1 + 2*eps)*s12^3*s13*(s12 + s13)*s23)']

#################################################################################################################################################################
###########################################################   FINAL COEFFICIENTS TO USE     ##################################################################### 
#################################################################################################################################################################

coeff = [ 
    
    # P1  - g1(mtsq)  coefficient
    [Coefficient(P1num,P1den,parameters),
    
    # P2  - g2(s12,mtsq) coefficient
    Coefficient(P2num,P2den,parameters),

    # P3  - g2(s13,mtsq) coefficient
    Coefficient(P3num,P3den,parameters),

    # P4  - g3(s23,mtsq) coefficient
    Coefficient(P4num,P4den,parameters),

    # P5  - g4(mhsq,mtsq) coefficient
    Coefficient(P5num,P5den,parameters),

    # P6  - g5(s12,mtsq) coefficient
    Coefficient(P6num,P6den,parameters),

    # P7  - g5(s13,mtsq) coefficient
    Coefficient(P7num,P7den,parameters),

    # P8  - g6(s23,mtsq) coefficient
    Coefficient(P8num,P8den,parameters),

    # P9  - g7(s12,mhsq,mtsq)  coefficient
    Coefficient(P9num,P9den,parameters),

    # P10  - g7(s13,mhsq,mtsq) coefficient
    Coefficient(P10num,P10den,parameters),

    # P11  - g8(s23,mhsq,mtsq) coefficient
    Coefficient(P11num,P11den,parameters),

    # P12  - g9(s12,s23,mhsq,mtsq) coefficient
    Coefficient(P12num,P12den,parameters),

    # P13  - g9(s12,s13,mhsq,mtsq) coefficient
    Coefficient(P13num,P13den,parameters),

    # P14  - g9(s23,s13,mhsq,mtsq) coefficient
    Coefficient(P14num,P14den,parameters)]

]

















