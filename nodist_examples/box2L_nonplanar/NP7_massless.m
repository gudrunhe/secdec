(* example is the non-planar massless 2-loop box with 7 propagators *)

momlist={k1,k2};

proplist={{0,{1,5}},{0,{2,6}},{0,{1,2}},{0,{3,5}},{0,{3,6}},{0,{4,6}},{0,{4,5}}};

numerator={1};

powerlist=Table[1,{i,Length[proplist]}];

Dim=4-2*eps;

prefactor=-Gamma[3+2*eps];

ExternalMomenta = {p1,p2,p3,p4};
externallegs=4;
KinematicInvariants = {s,t};
Masses={};

ScalarProductRules = {
  SP[p1,p1]->0,
  SP[p2,p2]->0,
  SP[p3,p3]->0,
  SP[p4,p4]->0,
  SP[p3,p2]->t/2,
  SP[p1,p3]->-t/2-s/2,
  SP[p1,p2]->s/2,
  SP[p1,p4]->t/2,SP[p2,p4]->-t/2-s/2,SP[p3,p4]->s/2};

 
