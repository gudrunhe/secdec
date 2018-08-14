(* case 8 of hep-ph/9605392, eq. (5.38) in Smirnov's book "Applied Asymptotic Expansions", expanded in (s/msq)^n up to n=2 *)

F55[s_,msq_,nmax_]:= Module[{prefac,delta,splusidelta,ll,deno,P2n0,P2n1,P2n2,P1n0,P1n1,P1n2,P0n0,P0n1,P0n2,res,rho,rr,rrexpanded},
                    
   prefac=-Exp[2*eps*EulerGamma];
   (* prefac is inverse of extracted prefac because pysd integral is without any prefac *)
   (* however the prefactor Gamma[2+2*eps] from the Feynman parametrisation is not included either in numerical result *)
   prefac=prefac/Gamma[2+2*eps];
   deno=s*msq;                 
   delta=10^(-15);  
   splusidelta=s+I*delta;                 
   ll=Log[-splusidelta/msq];
   P2n0 = -ll+1;
   P1n0 =  3/2*ll^2 - 3*ll - Pi^2/6 + 3;            
   P0n0 = -7/6*ll^3+7/2*ll^2+(Pi^2/3-7)*ll-2*Zeta[3]-1/3*Pi^2+7;
   P2n1 = -1/2*ll + 1/4; 
   P1n1 =  3/4*ll^2 - 9/4*ll - Pi^2/12 + 5/8;                   
   P0n1 = -7/12*ll^3+17/8*ll^2+(Pi^2/6-49/8)*ll-Zeta[3]-2/3*Pi^2+89/16;
   P2n2 = -1/3*ll + 1/9;                 
   P1n2 =  1/2*ll^2 - 11/6*ll - Pi^2/18 + 31/36;                   
   P0n2 = -7/18*ll^3+59/36*ll^2+(Pi^2/9-160/27)*ll-2/3*Zeta[3]-67/108*Pi^2+4709/648;
                    
   res[0]=P2n0/eps^2+P1n0/eps+P0n0;
   res[1]=P2n1/eps^2+P1n1/eps+P0n1;
   res[2]=P2n2/eps^2+P1n2/eps+P0n2;
   rho = -s/msq;                 
   rr=999;                 
   If[nmax<3,                 
      rr=prefac/deno*Sum[rho^n*res[n],{n,0,nmax}],
      Print["value for nmax must be in {0,1,2}"];
    ];
   rrexpanded=Normal[Series[rr,{eps,0,0}]];
   rrexpanded=N[rrexpanded];                     
   Return[rrexpanded]                   
]                    

(* full result for double pole *)

P2[s_,msq_]:= Module[{splusidelta},
splusidelta=s+I*10^(-15);
res=1/s^2 *(Log[-splusidelta/msq]*Log[1+splusidelta/msq]+PolyLog[2,-splusidelta/msq]);
Return[res]
]

test1=F55[1/2,1,0];
P2fulltest1=P2[0.5,1];

test2=F55[0.003,1,2];
P2fulltest2=P2[0.003,1];

ftprefac=Gamma[1-2*eps]/(Gamma[1+eps]^2*Gamma[1-eps]^2);
(* difference between ftprefac and Exp[2*eps*EulerGamma] is of order eps^3 *)


