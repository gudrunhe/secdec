(* funcs 43-46 and 91-94 *)

resultr43 = -0.006527348296
errorr43 = 6.053338122e-05

resultim43 = -1.105568539e-07
errorim43 = 6.10012437e-05


resultr44 = -0.004712835248
errorr44 = 2.096252991e-06

resultim44 = -1.049038956e-11
errorim44 = 8.46715336e-07


resultr45 = -0.0210702308
errorr45 = 2.094826794e-05

resultim45 = 2.747606374e-08
errorim45 = 1.043542297e-05

resultr46 = -0.001975309623
errorim46 = 1.777365208e-07

resultim46 = 1.731085609e-10
errorim46 = 1.256860875e-07

resultr91 = -0.007280845985
errorr91 = 7.132945549e-06

resultim91 = -7.763520058e-09
errorim91 = 5.451066669e-06

resultr92 = -0.0007090541307
errorr92 = 4.852243958e-07
resultim92 = -2.856948553e-12
errorim92 = 1.156976135e-07

resultr93 = -0.007280718368
errorr93 = 7.132530522e-06
resultim93 = 2.438998831e-09
errorim93 = 5.451109135e-06

resultr94 = -0.0006687134694
errorr94 = 4.484972986e-07
resultim94 = -7.021978039e-10
errorim94 = 2.073219368e-07

ReP4=resultr43+resultr44+resultr45+resultr46+resultr91+resultr92+resultr93+resultr94;
(* ReP4 = -0.0502250559201 *)

(* need to add results from 3l0h0 => eps^-4 when prefactoris taken into account *)

(* Stephen's analytic formula *)

myPolyLog[i__, x_?NumericQ] := Module[{},
   If[x > 0 && x < 1,
    PolyLog[i, x], Indeterminate]
   ];
myLog[x_] := Module[{},
   If[x > 0,
    Log[x], Indeterminate]
   ];
   
 b7NP[ins_, int_] := 
 Module[{FtRegion1, FtRegion2, FuRegion1, FuRegion2, T, U, V, inu},
  inu = -ins - int;
  (* Region i *)
  FtRegion1[s_, t_, u_] := ((
      (s^(-2*eps)) (-2/eps^4 + (37*Pi^4)/40 + (2*T^4)/3 + 96*U + 
         48*T*U - 12*T^2*U + (8*T^3*U)/3 + 12*T*U^2 - T^2*U^2 - 
         4*U^3 - (4*T*U^3)/3 + (4*U^4)/
          3 + (5*T + 7*U)/(2*eps^3) + ((-5*Pi^2)/12 - T^2 + 6*U - 
            4*T*U - U^2)/eps^2 + (10*Pi^2 - 48*T + 17*T^2 + 12*T*U)*
          myPolyLog[2, -(t/s)] + (48 - 60*T - 12*U)*
          myPolyLog[3, -(t/s)] + 
         86*myPolyLog[4, -(t/s)] + (-24 + 28*T - 6*U)*
          myPolyLog[1, 2, -(t/s)] - 26*myPolyLog[1, 3, -(t/s)] - 
         36*myPolyLog[2, 2, -(t/s)] + 
         I*Pi*(96 + 2/eps^3 + 48*T + 12*T^2 + (10*T^3)/3 - 24*T*U + 
            2*U^3 + (6 - T + U)/
             eps^2 - (24 + (31*Pi^2)/6 + 12*T + 2*T^2 + 2*T*U + 
               2*U^2)/eps + (-24 + 14*T + 18*U)*
             myPolyLog[2, -(t/s)] - 32*myPolyLog[3, -(t/s)] + 
            44*myPolyLog[1, 
              2, -(t/s)] + (Pi^2) (-2 + (28*T)/3 - U/3) - 
            89*Zeta[3]) + (Pi^2) ((-5*T^2)/6 - 
            14*U + (38*T*U)/3 + (25*U^2)/6) - (13*T + 33*U)*
          Zeta[3] + (-T^3/3 - 24*U - 12*T*U + 3*T*U^2 - U^3 - 
            2*T*myPolyLog[2, -(t/s)] + 2*myPolyLog[3, -(t/s)] - 
            2*myPolyLog[1, 2, -(t/s)] + (Pi^2) (5*T - 29*U)/
              6 + (19*Zeta[3])/2)/eps)
      ) /. {T -> myLog[-t/s], U -> myLog[-u/s], V -> myLog[-u/t]});
  FuRegion1[s_, t_, u_] := FtRegion1[s, u, t];
  (* Region ii *)
  FtRegion2[s_, t_, u_] := ((
      (t^(-2*eps)) (-2/eps^4 - (311*Pi^4)/120 + 96*T - 
         4*T^3 - (11*T^4)/6 + 96*V - 12*T^2*V - (13*T^3*V)/3 + 
         T^2*V^2 - 4*V^3 + 2*T*V^3 + 
           (4*V^4)/3 + (2*T + (7*V)/2)/
          eps^3 + ((31*Pi^2)/12 + 6*T + 2*T^2 + T*V + 6*V - V^2)/
          eps^2 + 
           (13*Pi^2 + 24*T + 6*T^2 - 18*T*V)*
          myPolyLog[2, -(s/t)] + (24 + 12*T - 18*V)*
          myPolyLog[3, -(s/t)] + 12*myPolyLog[4, -(s/t)] + 
           (24 - 44*T + 6*V)*myPolyLog[1, 2, -(s/t)] + 
         26*myPolyLog[1, 3, -(s/t)] - 62*myPolyLog[2, 2, -(s/t)] + 
           
         I*
          Pi (-5/(2*eps^3) - 48*T + 12*T^2 + (11*T^3)/3 - 48*V + 
            24*T*V + 2*T^2*V - 
            12*V^2 + (4*V^3)/3 + (T + 4*V)/eps^2 + 
                 ((5*Pi^2)/2 + 12*T + 4*T^2 + 12*V + 2*T*V - 3*V^2 - 
               2*myPolyLog[2, -(s/t)])/eps + (-48 + 14*T + 12*V)*
             myPolyLog[2, -(s/t)] + 
                 32*myPolyLog[3, -(s/t)] + 
            28*myPolyLog[1, 2, -(s/t)] + (Pi^2) (4 - 6*V) - 
            15*Zeta[3]) + 
           (Pi^2) (-14*T - (22*T^2)/3 - 
            14*V - (23*T*V)/3 + (37*V^2)/6) - (24 + 45*T + 39*V)*
          Zeta[3] + 
           (-24*T - (2*T^3)/3 - 24*V - 2*T^2*V - 2*T*V^2 - V^3 + 
            2*myPolyLog[1, 2, -(s/t)] - (Pi^2) (23*T + 41*V)/
              6 + (15*Zeta[3])/2)/eps)
      ) /. {T -> myLog[-t/s], U -> myLog[-u/s], V -> myLog[-u/t]});
  FuRegion2[s_, t_, u_] := ((
      (t^(-2*eps)) (-2/eps^4 - (311*Pi^4)/120 + 96*T - 
         4*T^3 - (11*T^4)/6 + 48*T*V - (5*T^3*V)/3 + 12*T*V^2 + 
         4*T^2*V^2 + (10*T*V^3)/3 + (2*V^4)/3 + 
           (2*T + (5*V)/2)/
          eps^3 + ((31*Pi^2)/12 + 6*T + 2*T^2 - T*V - V^2)/
          eps^2 + (11*Pi^2 - 24*T - 6*T^2 - 18*T*V)*
          myPolyLog[2, -(s/t)] - 
           (24 + 12*T + 18*V)*myPolyLog[3, -(s/t)] - 
         12*myPolyLog[4, -(s/t)] + (48 - 32*T + 26*V)*
          myPolyLog[1, 2, -(s/t)] + 
           86*myPolyLog[1, 3, -(s/t)] - 50*myPolyLog[2, 2, -(s/t)] + 
           
         I*Pi (-96 - 7/(2*eps^3) + (7*T^3)/3 - 48*V - 4*T^2*V - 
            12*V^2 - 2*T*V^2 - (4*V^3)/3 + (-6 - T + 4*V)/eps^2 + 
                 (24 + (11*Pi^2)/2 + 2*T^2 + 12*V + 2*T*V - V^2 - 
               2*myPolyLog[2, -(s/t)])/eps + (-24 + 26*T - 8*V)*
             myPolyLog[2, -(s/t)] + 
                 44*myPolyLog[3, -(s/t)] - 
            24*myPolyLog[1, 2, -(s/t)] + (Pi^2) (6 + 6*T - 8*V) + 
            75*Zeta[3]) + 
           (Pi^2) (-14*T - (22*T^2)/3 - (7*T*V)/3 + (31*V^2)/
             6) - (24 + 45*T + 21*V)*Zeta[3] + 
           (-24*T - (2*T^3)/3 - 12*T*V - 4*T^2*V - 2*T*V^2 - V^3/3 + 
            2*myPolyLog[1, 2, -(s/t)] - (Pi^2) (23*T + 31*V)/
              6 + (15*Zeta[3])/2)/eps)
      ) /. {T -> myLog[-t/s], U -> myLog[-u/s], V -> myLog[-u/t]});
  Series[Piecewise[
    {
     {Gamma[1 + eps]^2*(FtRegion1[ins, int, inu]/(ins^2*int) + 
         FuRegion1[ins, int, inu]/(ins^2*inu)), (inu < 0 && int < 0 &&
         ins > 0)},
     {Gamma[1 + eps]^2*(FtRegion2[ins, int, inu]/(ins^2*int) + 
         FuRegion2[ins, int, inu]/(ins^2*inu)), (inu < 0 && ins < 0 &&
         int > 0)},
     {Gamma[1 + eps]^2*(FtRegion2[ins, inu, int]/(ins^2*inu) + 
         FuRegion2[ins, inu, int]/(ins^2*int)), (ins < 0 && int < 0 &&
         inu > 0)} (* t <> u*)
     }, Indeterminate
    ], {eps, 0, 0}]
  ];  
  
 nanalytic= N[Normal[b7NP[9, -2.5]]];
  
(*
nanalytic = 
     0.013675213675213675/eps^4 - (0.04438309817390562 + 
       0.04296195081832196*I)/eps^3 + 
     (0.1339431853817464 + 0.10038251452523675*I)/eps^2 - 
     (0.6971699092107224 - 1.475821473466632*I)/eps +
     (-0.050022596625423854 - 4.271109175940437*I)  
*)
