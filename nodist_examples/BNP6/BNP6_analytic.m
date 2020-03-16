(* from hep-ph/9909506 *)

BNP6[s_,t_,u_]:=Module[{prefac,T,U,part1,part2,res},
prefac=3/eps^2*Gamma[1+2*eps]*Gamma[1-eps]^3/Gamma[1-3*eps]/(1+4*eps);
T=Log[-t/s];
U=Log[-u/s];
part1=s^(-2*eps)/(s*u)*((T+I*Pi)*(1-2*eps*U)+
      2*eps^2*(2*PolyLog[3,-t/s]-4*PolyLog[1,2,-t/s]+2*Zeta[3]-T^3/3+
      T*(-2*PolyLog[2,-t/s]+U^2-Pi^2)+I*Pi*(2*PolyLog[2,-t/s]+U^2-Pi^2/3) ));
part2=s^(-2*eps)/(s*t)*((U+I*Pi)*(1-2*eps*T)+
      2*eps^2*(2*PolyLog[3,-u/s]-4*PolyLog[1,2,-u/s]+2*Zeta[3]-U^3/3+
      U*(-2*PolyLog[2,-u/s]+T^2-Pi^2)+I*Pi*(2*PolyLog[2,-u/s]+T^2-Pi^2/3) ));
res=prefac*(part1+part2);
Return[res];
];

test1=N[Normal[Series[BNP6[9,-2.5,-6.5], {eps,0,0}]]];
testEuc=N[Normal[Series[BNP6[-9,-2.5,-6.5], {eps,0,0}]]];

(*
test1 = (0.7966721383373115 - 18.22048104236002*I) + 
     (0.10907856854318447 - 0.5799863360473464*I)/eps^2 - 
     (0.8876663743916553 - 4.360251717854891*I)/eps

testEuc = (-105.05122507795315 - 74.2320300324536*I) - 
     (0.10907856854318447 - 1.1599726720946928*I)/eps^2 + 
     (15.464312874462092 - 7.3497817173123465*I)/eps
*)

