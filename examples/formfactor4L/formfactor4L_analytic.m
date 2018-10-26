(* result from equation (7.1) in arXiv:1510.06758 *)

zeta53 = NSum[Sum[n^-3 m^-5, {n, 1, m - 1}], {m, 1, \[Infinity]}, AccuracyGoal -> 30, WorkingPrecision -> 100];

formfactor4L = N[
 18/5 * Zeta[2]^2  Zeta[3] -
  5 Zeta[2] Zeta[
    5] + (24 Zeta[2] Zeta[3] + 20 Zeta[5] - 188/105 Zeta[2]^3 -
     17 Zeta[3]^2 + 9 Zeta[2]^2  Zeta[3] - 47 Zeta[2] Zeta[5] -
     21 Zeta[7] + 6883/2100 *Zeta[2]^4 + 49/2*Zeta[2] Zeta[3]^2 +
     1/2 Zeta[3] Zeta[5] - 9 zeta53)*eps , 20]

(* formfactor4L = 3.1807380843134699650 + 46.104303262308462367 eps *)
