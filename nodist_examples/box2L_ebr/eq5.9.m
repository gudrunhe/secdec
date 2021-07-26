(* ::Package:: *)

<<HPL`


integral[ep0_,s0_,t0_,u0_,mt0_]:=
Module[{ep=ep0,s=s0,t=t0,mt1=mt0,mus2,v,mt,im},
	v=-t/s;
	mt=mt1/Sqrt[s];
	im=I;
	mus2=Exp[2 im Pi ep] Exp[-2 ep EulerGamma]/ s^(3+2 ep);
	mus2*((4*im*Pi*(12*HPL[{3}, v] - 6*HPL[{2}, v]*Log[v] + 
     (Log[mt^2] - Log[v])*(Pi^2 + Log[mt^2]^2 - 2*Log[mt^2]*Log[v] + 
       Log[v]^2)))/(3*v) - (7*Pi^4 + 360*HPL[{4}, v] + 10*Pi^2*Log[mt^2]^2 - 
    15*Log[mt^2]^4 - 240*HPL[{3}, v]*Log[v] - 40*Pi^2*Log[mt^2]*Log[v] + 
    40*Log[mt^2]^3*Log[v] + 30*Pi^2*Log[v]^2 + 60*HPL[{2}, v]*Log[v]^2 - 
    30*Log[mt^2]^2*Log[v]^2 + 5*Log[v]^4 + 60*Log[mt^2]*Zeta[3] - 
    60*Log[v]*Zeta[3])/(15*v) + 
  ep*((im*Pi*(61*Pi^4 - 360*HPL[{2}, v]^2 - 1440*HPL[{4}, v] - 
       1080*HPL[{2, 2}, v] + 1080*HPL[{3}, v]*Log[mt^2] - 
       30*Pi^2*Log[mt^2]^2 - 315*Log[mt^2]^4 + 2160*HPL[{3}, v]*Log[1 - v] - 
       2880*HPL[{3}, v]*Log[v] + 720*HPL[{2, 1}, v]*Log[v] - 
       420*Pi^2*Log[mt^2]*Log[v] + 780*Log[mt^2]^3*Log[v] - 
       180*Pi^2*Log[1 - v]*Log[v] - 180*Log[mt^2]^2*Log[1 - v]*Log[v] + 
       450*Pi^2*Log[v]^2 - 450*Log[mt^2]^2*Log[v]^2 + 
       180*Log[mt^2]*Log[1 - v]*Log[v]^2 - 180*Log[mt^2]*Log[v]^3 - 
       60*Log[1 - v]*Log[v]^3 + 165*Log[v]^4 + 60*HPL[{2}, v]*
        (5*Pi^2 - 3*Log[mt^2]^2 - 6*Log[mt^2]*Log[v] - 18*Log[1 - v]*Log[v] + 
         27*Log[v]^2) - 2160*Log[1 - v]*Zeta[3]))/(90*v) + 
    (7560*HPL[{5}, v] + 1440*HPL[{2, 3}, v] + 2520*HPL[{3, 2}, v] + 
      218*Pi^4*Log[mt^2] + 120*Pi^2*HPL[{2}, v]*Log[mt^2] - 
      1440*HPL[{4}, v]*Log[mt^2] + 200*Pi^2*Log[mt^2]^3 - 120*Log[mt^2]^5 + 
      36*Pi^4*Log[1 - v] - 3240*HPL[{4}, v]*Log[1 - v] - 157*Pi^4*Log[v] - 
      540*Pi^2*HPL[{2}, v]*Log[v] - 360*HPL[{2}, v]^2*Log[v] - 
      1080*HPL[{2, 2}, v]*Log[v] - 630*Pi^2*Log[mt^2]^2*Log[v] - 
      180*HPL[{2}, v]*Log[mt^2]^2*Log[v] + 285*Log[mt^2]^4*Log[v] + 
      120*Pi^2*Log[mt^2]*Log[1 - v]*Log[v] + 360*HPL[{2, 1}, v]*Log[v]^2 + 
      420*Pi^2*Log[mt^2]*Log[v]^2 - 180*Log[mt^2]^3*Log[v]^2 - 
      150*Pi^2*Log[1 - v]*Log[v]^2 - 540*HPL[{2}, v]*Log[1 - v]*Log[v]^2 - 
      90*Log[mt^2]^2*Log[1 - v]*Log[v]^2 + 10*Pi^2*Log[v]^3 + 
      540*HPL[{2}, v]*Log[v]^3 + 30*Log[mt^2]^2*Log[v]^3 + 
      120*Log[mt^2]*Log[1 - v]*Log[v]^3 - 60*Log[mt^2]*Log[v]^4 - 
      45*Log[1 - v]*Log[v]^4 + 45*Log[v]^5 + 60*HPL[{3}, v]*
       (13*Pi^2 + 6*HPL[{2}, v] + 3*Log[mt^2]^2 + 12*Log[mt^2]*Log[v] + 
        36*Log[1 - v]*Log[v] - 30*Log[v]^2) + 540*Pi^2*Zeta[3] - 
      5040*HPL[{2}, v]*Zeta[3] - 360*Log[mt^2]^2*Zeta[3] + 
      720*Log[mt^2]*Log[v]*Zeta[3] - 2160*Log[1 - v]*Log[v]*Zeta[3] - 
      360*Log[v]^2*Zeta[3] - 4680*Zeta[5])/(90*v)))
]
s=4.0;
t=-2.82842712475;
u=-2.82842712475;
mt=Sqrt[0.1];
Print["Series expansion:"]
Series[integral[ep,s,t,u,mt],{ep,0,1}]//N
