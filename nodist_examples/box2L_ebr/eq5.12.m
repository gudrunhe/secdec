(* ::Package:: *)

<<HPL`


integral[ep0_,s0_,t0_,u0_,mt0_]:=
Module[{ep=ep0,s=s0,t=t0,mt1=mt0,mus2,v,mt,im},
	v=-t/s;
	mt=mt1/Sqrt[s];
	im=I;
	mus2=Exp[2 im Pi ep] Exp[-2 ep EulerGamma]/ (s)^(3+2 ep);
	mus2*(((2*im*Pi*(Log[mt^2] - Log[v]))/v - 
    (Pi^2 - 2*Log[mt^2]^2 + 2*Log[mt^2]*Log[v])/v)/ep^2 + 
  ((im*Pi*(2*Pi^2 - 2*HPL[{2}, v] - 5*Log[mt^2]^2 + 2*Log[mt^2]*Log[v] - 
       2*Log[1 - v]*Log[v] + 3*Log[v]^2))/v + 
    (6*HPL[{3}, v] + 10*Pi^2*Log[mt^2] - 4*Log[mt^2]^3 - 4*Pi^2*Log[v] - 
      6*HPL[{2}, v]*Log[v] - 3*Log[mt^2]^2*Log[v] + 6*Log[mt^2]*Log[v]^2 - 
      3*Log[1 - v]*Log[v]^2 + Log[v]^3 - 42*Zeta[3])/(3*v))/ep - 
  (im*Pi*(-60*HPL[{3}, v] + 6*HPL[{2, 1}, v] + 5*Pi^2*Log[mt^2] - 
     10*Log[mt^2]^3 - Pi^2*Log[1 - v] + 6*HPL[{2}, v]*Log[1 - v] + 
     7*Pi^2*Log[v] + 18*HPL[{2}, v]*Log[v] - 6*Log[mt^2]^2*Log[v] + 
     3*Log[1 - v]^2*Log[v] + 6*Log[mt^2]*Log[v]^2 - 12*Log[1 - v]*Log[v]^2 + 
     10*Log[v]^3 - 72*Zeta[3]))/(3*v) + 
  (259*Pi^4 + 90*HPL[{2}, v]^2 - 6120*HPL[{4}, v] - 180*HPL[{2, 2}, v] - 
    1050*Pi^2*Log[mt^2]^2 + 90*Log[mt^2]^4 + 360*HPL[{3}, v]*Log[1 - v] + 
    3600*HPL[{3}, v]*Log[v] - 360*HPL[{2, 1}, v]*Log[v] + 
    360*Pi^2*Log[mt^2]*Log[v] + 240*Log[mt^2]^3*Log[v] - 
    420*Pi^2*Log[1 - v]*Log[v] + 330*Pi^2*Log[v]^2 - 
    180*Log[mt^2]^2*Log[v]^2 - 90*Log[1 - v]^2*Log[v]^2 + 
    240*Log[1 - v]*Log[v]^3 - 150*Log[v]^4 - 60*HPL[{2}, v]*
     (7*Pi^2 + 6*Log[1 - v]*Log[v] + 9*Log[v]^2) + 1800*Log[mt^2]*Zeta[3] - 
    360*Log[1 - v]*Zeta[3] + 2520*Log[v]*Zeta[3])/(180*v) + 
  ep*(-(im*Pi*(94*Pi^4 + 540*HPL[{2}, v]^2 + 1440*HPL[{4}, v] + 
        2520*HPL[{2, 2}, v] + 360*HPL[{2, 1, 1}, v] - 1080*HPL[{3}, v]*
         Log[mt^2] - 600*Pi^2*Log[mt^2]^2 + 255*Log[mt^2]^4 - 
        5040*HPL[{3}, v]*Log[1 - v] + 360*HPL[{2, 1}, v]*Log[1 - v] - 
        30*Pi^2*Log[1 - v]^2 + 7200*HPL[{3}, v]*Log[v] - 
        1800*HPL[{2, 1}, v]*Log[v] + 480*Pi^2*Log[mt^2]*Log[v] + 
        420*Log[mt^2]^3*Log[v] + 240*Pi^2*Log[1 - v]*Log[v] + 
        180*Log[mt^2]^2*Log[1 - v]*Log[v] + 60*Log[1 - v]^3*Log[v] - 
        600*Pi^2*Log[v]^2 - 270*Log[mt^2]^2*Log[v]^2 - 
        180*Log[mt^2]*Log[1 - v]*Log[v]^2 - 360*Log[1 - v]^2*Log[v]^2 + 
        180*Log[mt^2]*Log[v]^3 + 780*Log[1 - v]*Log[v]^3 - 585*Log[v]^4 - 
        180*HPL[{2}, v]*(5*Pi^2 - Log[mt^2]^2 - Log[1 - v]^2 - 
          2*Log[mt^2]*Log[v] - 10*Log[1 - v]*Log[v] + 17*Log[v]^2) + 
        5280*Log[mt^2]*Zeta[3] + 4680*Log[1 - v]*Zeta[3] + 
        1920*Log[v]*Zeta[3]))/(180*v) - 
    (-17640*HPL[{5}, v] + 420*Pi^2*HPL[{2, 1}, v] - 
      60*HPL[{2}, v]*HPL[{2, 1}, v] - 3180*HPL[{2, 3}, v] - 
      5940*HPL[{3, 2}, v] + 60*HPL[{2, 1, 2}, v] + 180*HPL[{2, 2, 1}, v] - 
      25*Pi^4*Log[mt^2] - 120*Pi^2*HPL[{2}, v]*Log[mt^2] + 
      1440*HPL[{4}, v]*Log[mt^2] - 680*Pi^2*Log[mt^2]^3 + 24*Log[mt^2]^5 - 
      161*Pi^4*Log[1 - v] + 420*Pi^2*HPL[{2}, v]*Log[1 - v] - 
      90*HPL[{2}, v]^2*Log[1 - v] + 8280*HPL[{4}, v]*Log[1 - v] + 
      180*HPL[{2, 2}, v]*Log[1 - v] + 599*Pi^4*Log[v] + 
      900*Pi^2*HPL[{2}, v]*Log[v] + 540*HPL[{2}, v]^2*Log[v] + 
      2520*HPL[{2, 2}, v]*Log[v] + 360*HPL[{2, 1, 1}, v]*Log[v] - 
      240*Pi^2*Log[mt^2]^2*Log[v] + 180*HPL[{2}, v]*Log[mt^2]^2*Log[v] + 
      135*Log[mt^2]^4*Log[v] + 360*HPL[{2, 1}, v]*Log[1 - v]*Log[v] - 
      120*Pi^2*Log[mt^2]*Log[1 - v]*Log[v] + 210*Pi^2*Log[1 - v]^2*Log[v] + 
      180*HPL[{2}, v]*Log[1 - v]^2*Log[v] - 900*HPL[{2, 1}, v]*Log[v]^2 + 
      240*Pi^2*Log[mt^2]*Log[v]^2 - 60*Log[mt^2]^3*Log[v]^2 - 
      540*Pi^2*Log[1 - v]*Log[v]^2 + 900*HPL[{2}, v]*Log[1 - v]*Log[v]^2 + 
      90*Log[mt^2]^2*Log[1 - v]*Log[v]^2 + 30*Log[1 - v]^3*Log[v]^2 + 
      440*Pi^2*Log[v]^3 - 1020*HPL[{2}, v]*Log[v]^3 - 
      30*Log[mt^2]^2*Log[v]^3 - 120*Log[mt^2]*Log[1 - v]*Log[v]^3 - 
      120*Log[1 - v]^2*Log[v]^3 + 60*Log[mt^2]*Log[v]^4 + 
      225*Log[1 - v]*Log[v]^4 - 129*Log[v]^5 - 60*HPL[{3}, v]*
       (7*HPL[{2}, v] + 3*(16*Pi^2 + Log[mt^2]^2 + Log[1 - v]^2 + 
          4*Log[mt^2]*Log[v] + 28*Log[1 - v]*Log[v] - 22*Log[v]^2)) - 
      5220*Pi^2*Zeta[3] + 11520*HPL[{2}, v]*Zeta[3] + 
      3300*Log[mt^2]^2*Zeta[3] + 180*Log[1 - v]^2*Zeta[3] - 
      1320*Log[mt^2]*Log[v]*Zeta[3] + 4680*Log[1 - v]*Log[v]*Zeta[3] + 
      1620*Log[v]^2*Zeta[3] + 21780*Zeta[5])/(180*v)))
]

s=5.3;
t=-1.86;
u=0;
mt=Sqrt[0.1];
Print["Series expansion: "]
Series[integral[ep,s,t,u,mt],{ep,0,0}]//N



