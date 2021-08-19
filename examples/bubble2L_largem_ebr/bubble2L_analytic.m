(* Analytic result as expansion. Converges only for s << msq *)

I3[s_,msq_,nmax_]:= Module[{splusidelta,f0,f1},
splusidelta=s+I*10^(-15);
f0=3/2/msq*Sum[Factorial[n]*Gamma[3/2]/(Gamma[n+3/2]*(n+1)^2)*(-s/4/msq)^n,{n,0,nmax}];
f1=-1/2/msq*Sum[Factorial[n]*Gamma[3/2]/(Gamma[n+3/2]*(n+1))*(-s/4/msq)^n,{n,0,nmax}];
res=f0+f1*Log[-splusidelta/msq];
Return[res]
]

test=I3[0.002,4,0]
