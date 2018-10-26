(* result from equation (A3) in arXiv:1705.06483 *)

zeta35 = NSum[Sum[n^-3 m^-5, {n, 1, m - 1}], {m, 1, \[Infinity]}, AccuracyGoal -> 30, WorkingPrecision -> 100];

epsm2 = 147/16*Zeta[7];
epsm1 = -(147/16*Zeta[7]+27/2*Zeta[3]*Zeta[5]+27/10*zeta35-2063/504000*Pi^8);
eps0 = unknown;
zeta35 = NSum[Sum[n^-3 m^-5, {n, 1, m - 1}], {m, 1, \[Infinity]}, AccuracyGoal -> 30, WorkingPrecision -> 100]; (* approx. 0.037707673 *)
prefac = Gamma[1-eps]^2*Gamma[1+eps]/Gamma[2-2*eps];

(* pysecdec prefactor (from Feynman parametrisation) is Gamma[6*eps] *)
(* prefactor in 1705.06483 eq.(A3) is prefac=Gamma[1-eps]^2*Gamma[1+eps]/Gamma[2-2*eps] per loop *)

bubble6L = N[Normal[Series[prefac^6*(epsm2/eps^2+epsm1/eps), {eps,0,-1}]]]

(* bubble6L = 9.264208985946416/eps^2 + 91.73175282208716/eps *)
