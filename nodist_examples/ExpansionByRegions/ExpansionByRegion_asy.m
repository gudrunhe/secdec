(* ::Package:: *)

(* ::Input:: *)
(**)


(*
The path to the executable qhull (on your computer) can be set in the Options[QHull] in asy2.m.
*)
Get["./asy2.m"];

(*The variable "x" is used below as the collective name for the parametric-integration-variables, 
for it is being hard coded in ASY where UNFORTUNATELY it is also used as the expansion-parameter! 
This gives a bit trouble when "applying the region expansions" for purely technical reasons. 
For the time being, this is simply (albeit conceptually not elegantly) circumvented by an inner 
"renaming" of the expansion-parameter when applying "applying the region expansions". 
It is also easy to just modify the ASM-code to remove this unpleasant limitation.*)

ClearAll[ExpansionByRegion,ExpansionByRegionLO,SymanzikParaIntegrand,x,Dr,secdec];

Protect[x,Dr];
(*
The exkinepowerscaling is assumed to be of the form of a list of purely monic powers of the unique scaling(expansion)-parameter.

The last non-optional argument "po" specifies the (absolute) power-order in "x" to which the loop-integral will be Taylor-expanded.
The absoute power-orders of the "relative leading-power-terms" are listed in the first entry of the returned list, 
while the second entry gives the list of the expanded/truncated loop-integrands per region. 
*)

ExpansionByRegion[Loopmomenta_List,FPropagators_List,FPdenPowers_List,kinereplacement_List,exkinepowerscalingraw_List,po_Integer,OptionsPattern[{ExpansionSpaceDim-> 4,KeepEPRarameterAs-> 1,AnalyticRegulator-> apr,UFlist-> False,Verbosity-> True}]]:=Module[{svx,sv,Sup,Sfp,nls,lcmN,exkinepowerscaling,Xregpowerraw,Xregpower,Xregpowerscaling,SupXreg,SfpXreg,LOpowerXregion,sign,gammaPrefactor,deltafactor,xweightfactor,redFUintegrandList,redFUintegrandListTaylorExpanded,SymanzikParaIntegrandList},
(***************************A few preliminary check of input arguments***********************************************************)
With[{svmonolist=Map[Last,exkinepowerscalingraw]},
If[(Union[(svmonolist/.{x-> 1})]=!={1}),
Print["The power-scaling behaviors of all external-kinematics should be specified such that the right-hand sides of each rule
 is purely a monic power of the unique scaling-parameter."];Abort[];]];
(*identify the scaling(expansion)-parameter from the input "exkinepowerscalingraw".*)
svx=With[{svfind=Complement[Variables[exkinepowerscalingraw/.{Rule-> List}],Variables[First[Transpose[exkinepowerscalingraw/.{Rule-> List}]]]]},
If[Length[svfind]===1,First[svfind],Print["The input specifying the power-scaling behavior of external-kinematics is not given in the proper format, and the scaling parameter is not found."];Abort[]]];
If[Not[SubsetQ[{OptionValue[AnalyticRegulator]},Variables[FPdenPowers]]],Print["The list of propagator powers contain symbols other than the specified analytic-regulator (its default name is \"apr\" set in the option 
\"AnalyticRegulator\")!"];Abort[]];
If[OptionValue[Verbosity]===True,Print["The scaling-parameter appearing in the input is : ",svx]];
(**************************Get the region-vectors from ASY*****************************************************)
(*get the UF-polynomials*)
{Sup,Sfp,nls}=UF[Loopmomenta,FPropagators,kinereplacement];
(*get the list of power-scaling vectors for all scaleful expansion-regions*)
Xregpowerraw=AlphaRepExpand[Loopmomenta,FPropagators,kinereplacement,exkinepowerscalingraw,PreResolve->False];
If[Not[FreeQ[Xregpowerraw,{}]],Print["No region is found, and the expansion is aborted."];Abort[]];
(*Obtain the integer which is needed to make sure all power-scaling vectors have only integer components, 
which in most considered cases is just 1. Multiply this lcmN to all region-vectors, including the one for external-kinematics
amounts to a global change in the definition of the power-scaling parameter "sv".*)
lcmN=With[{denfs=DeleteCases[Abs[Flatten[Xregpowerraw]],0]},If[denfs==={},1,Apply[LCM,Map[Denominator,denfs]]]];
If[OptionValue[Verbosity]===True,Print["The least-common-multipler of all denomaintors of region-vectors' components is: ",lcmN]];
(*for the time-being, just rename the expansion-parameter "x" (which is saved in svx) as "sv"*)
exkinepowerscaling=exkinepowerscalingraw/.{Rule[l_,r_]:>Rule[l,r^lcmN] }/.{svx-> sv};
If[OptionValue[Verbosity]===True,Print["The power-scaling of external-kinematics (in dummay x) used now reads as: ",exkinepowerscaling/.{sv-> x}]];
(**************************Apply the expansion-by-regions (Begin)*****************************************************)
(*a uniform shift such that the minimum scaling power is 0, locally within each expansion-region*)
Xregpower=Map[(#-Min[#])&,(lcmN*Xregpowerraw)];
If[OptionValue[Verbosity]===True,Print["The list of (reformulated) region-vectors reads as: ",Xregpower]];
Xregpowerscaling=Map[Thread[Rule[Table[x[i],{i,1,Length[#]}],Power[sv,#]]]&,Xregpower];
(*************************************)
(*Manifest the power-scaling behavior of each term of UF-polynomials in each expansion-region (no expansion or truncation yet).
We pull out a pure monomial power factor of the expansion-parameter such that the thus-reduced polynomial has a "constant" term.
The formal pure monomial power factor pulled out could be sv^0 = 1, which is fine.*)
SupXreg=With[{SupTerms=MonomialList[Sup],SupMonomials=(MonomialList[Sup]/.{Times[nc_Integer,rest__]:> Times[1,rest]})},Map[(With[{svPSlist=Map[(Exponent[#,sv,Min])&,SupMonomials/.exkinepowerscaling/.##]},{Min[svPSlist],Total[SupTerms*Power[sv,(svPSlist-Min[svPSlist])]]}])&,Xregpowerscaling]];
SfpXreg=With[{SfpTerms=MonomialList[Sfp],SfpMonomials=(MonomialList[Sfp]/.{Times[nc_Integer,rest__]:> Times[1,rest]})},Map[(With[{svPSlist=Map[(Exponent[#,sv,Min])&,SfpMonomials/.exkinepowerscaling/.##]},{Min[svPSlist],Total[SfpTerms*Power[sv,(svPSlist-Min[svPSlist])]]}])&,Xregpowerscaling]];
(*************************************)
LOpowerXregion=MapIndexed[(Dot[FPdenPowers,#1]+(First[SupXreg[[First[#2]]]]*(Total[FPdenPowers]-(nls+1)*Dr/2))-(First[SfpXreg[[First[#2]]]]*(Total[FPdenPowers]-(nls)*Dr/2)))&,Xregpower]/.{Dr-> OptionValue[ExpansionSpaceDim],OptionValue[AnalyticRegulator]-> 0}//Simplify;
If[OptionValue[Verbosity]===True,Print["The LO-power per region: ",LOpowerXregion]];
(*************************************)
sign=(-1)^(Total[FPdenPowers]);
gammaPrefactor=Gamma[Total[FPdenPowers]-nls*Dr/2]/Product[Gamma[FPdenPowers[[i]]],{i,1,Length[FPdenPowers]}];
deltafactor=DiracDelta[1-Sum[x[i],{i,1,Length[FPropagators]}]];
xweightfactor=Product[Power[x[i],FPdenPowers[[i]]-1],{i,1,Length[FPdenPowers]}];
(*compose the core piece of the sv-manifested FU-parametric loop-integrands of resulting loop-integrals in each region.
The only piece suppressed is just the common prefactor "gammaPrefactor*deltafactor*xweightfactor".*)
redFUintegrandList=Table[((Last[SupXreg[[n]]])^(Total[FPdenPowers]-(nls+1)*Dr/2))*(Last[SfpXreg[[n]]])^(-(Total[FPdenPowers]-nls*Dr/2)),{n,1,Length[Xregpowerraw]}];
(*The (po*lcmN-LOpowerXregion[[First[#2]]]) can be negative for some regions, and then Series[] gives 0 as expected.*)
redFUintegrandListTaylorExpanded=MapIndexed[CompoundExpression[If[OptionValue[Verbosity]===True,Print["The Taylor-expansion-depth w.r.t the LO-result in region-",First[#2]," is ",(po*lcmN-LOpowerXregion[[First[#2]]])+1]],Normal[Series[#1,{sv,0,(po*lcmN-LOpowerXregion[[First[#2]]])}]]]&,redFUintegrandList];
(*The following definition agrees with SecDec convention*)
SymanzikParaIntegrandList=sign*gammaPrefactor*deltafactor*xweightfactor*((Power[sv,LOpowerXregion]*redFUintegrandListTaylorExpanded)//.{sv-> OptionValue[KeepEPRarameterAs]});
(**************************Apply the expansion-by-regions (End)*****************************************************)
Return[If[OptionValue[UFlist],{SupXreg,SfpXreg},{(LOpowerXregion/lcmN),SymanzikParaIntegrandList}]];
];


(*
Strictly speaking this function keeps the "relative leading-power-term" from every non-vanishing expansion-region, while these terms 
may not be of the same absoute power-order. The absoute power-orders of these "relative leading-power-terms" are listed in the 
first entry of the returned list, while the second entry gives the list of the expanded/truncated loop-integrands per region. 
*)

ExpansionByRegionLO[Loopmomenta_List,FPropagators_List,FPdenPowers_List,kinereplacement_List,exkinepowerscaling_List,OptionsPattern[{UFlist-> False}]]:=Module[{sv,Sup,Sfp,nls,Xregionalpowerraw,Xregionalpower,Xregionalpowerscaling,SupXregLE,SfpXregLE,LOpowerXregion,sign,gammaPrefactor,deltafactor,xweightfactor,SymanzikParaIntegrandList},
(*identify the scaling(expansion)-parameter from the input "exkinepowerscaling".*)
sv=With[{svfind=Complement[Variables[exkinepowerscaling/.{Rule-> List}],Variables[First[Transpose[exkinepowerscaling/.{Rule-> List}]]]]},
If[Length[svfind]===1,First[svfind],Print["The input specifying the power-scaling behavior of external-kinematics is not given in the proper format, and the scaling parameter is not found."];Abort[]]];
(*get the UF-polynomials*)
{Sup,Sfp,nls}=UF[Loopmomenta,FPropagators,kinereplacement];
(*get the list of power-scaling vectors for all scaleful expansion-regions*)
Xregionalpowerraw=AlphaRepExpand[Loopmomenta,FPropagators,kinereplacement,exkinepowerscaling,PreResolve->False];
If[Not[FreeQ[Xregionalpowerraw,{}]],Print["No region is found, and the expansion is aborted."];Abort[]];
(*a uniform shift such that the minimum scaling power is 0*)
Xregionalpower=Map[(#-Min[#])&,Xregionalpowerraw];
Xregionalpowerscaling=Map[Thread[Rule[Table[x[i],{i,1,Length[#]}],Power[sv,#]]]&,Xregionalpower];
(*get the Leading-Power-Order terms in the Taylor-expanded UF-polynomials in each expansion region*)
SupXregLE=With[{SupTerms=MonomialList[Sup]},Map[(With[{svPSlist=Map[(Exponent[#,sv])&,SupTerms/.exkinepowerscaling/.#]},{Min[svPSlist],Total[Pick[SupTerms,svPSlist,Min[svPSlist]]]}])&,Xregionalpowerscaling]];
SfpXregLE=With[{SfpTerms=MonomialList[Sfp]},Map[(With[{svPSlist=Map[(Exponent[#,sv])&,SfpTerms/.exkinepowerscaling/.#]},{Min[svPSlist],Total[Pick[SfpTerms,svPSlist,Min[svPSlist]]]}])&,Xregionalpowerscaling]];
(****************************************************************************************)
LOpowerXregion=MapIndexed[(Dot[FPdenPowers,#1]+(First[SupXregLE[[First[#2]]]]*(Total[FPdenPowers]-(nls+1)*Dr/2))-(First[SfpXregLE[[First[#2]]]]*(Total[FPdenPowers]-(nls)*Dr/2)))&,Xregionalpower]//Simplify;
(*compose the FU-parametric loop-integrands of resulting loop-integrals in each region, according to SecDec convention*)
sign=(-1)^(Total[FPdenPowers]);
gammaPrefactor=Gamma[Total[FPdenPowers]-nls*Dr/2]/Product[Gamma[FPdenPowers[[i]]],{i,1,Length[FPdenPowers]}];
deltafactor=DiracDelta[1-Sum[x[i],{i,1,Length[FPropagators]}]];
xweightfactor=Product[Power[x[i],FPdenPowers[[i]]-1],{i,1,Length[FPdenPowers]}];
SymanzikParaIntegrandList=(sign*gammaPrefactor*deltafactor*xweightfactor)*(Table[((Last[SupXregLE[[n]]])^(Total[FPdenPowers]-(nls+1)*Dr/2))*(Last[SfpXregLE[[n]]])^(-(Total[FPdenPowers]-nls*Dr/2)),{n,1,Length[Xregionalpowerraw]}]);
(****************************************************************************************)
Return[If[OptionValue[UFlist],{SupXregLE,SfpXregLE},{LOpowerXregion/.{Dr-> 4},SymanzikParaIntegrandList}]];
];

SymanzikParaIntegrand[LoopList_,PropagatorList_,PropagatorPowerList_,KeniRules_,paraV_: x,stdim_: Dr,FIConvention_: secdec]:=Module[{totPP,sign,pifactor,gammaPrefactor,deltafactor,xweightfactor,Ufactor,Ffactor,loopNum,propNum,FIntegrad},
{Ufactor,Ffactor,loopNum}=UF[LoopList,PropagatorList,KeniRules];
propNum=Length[PropagatorList];
totPP=Total[PropagatorPowerList];
pifactor=(Pi)^(loopNum*stdim/2);
sign=(-1)^(totPP);
gammaPrefactor=Gamma[totPP-loopNum*stdim/2]/Product[Gamma[PropagatorPowerList[[i]]],{i,1,propNum}];
deltafactor=DiracDelta[1-Sum[paraV[i],{i,1,propNum}]];
xweightfactor=Product[Power[paraV[i],PropagatorPowerList[[i]]-1],{i,1,propNum}];
FIntegrad=sign*gammaPrefactor*deltafactor*xweightfactor*(Ufactor^(totPP-(loopNum+1)*stdim/2))*Ffactor^(-(totPP-loopNum*stdim/2));
Return[If[FIConvention===secdec,FIntegrad,pifactor*FIntegrad]];
];



ClearAll[ExpansionByRegionGDRI];

(*
The possible non-integer symbol appearing in the power-exponents of the polynomials of the integrand is restricted to be named as "Dr", which is the spacetime-dimension for DR-loop-integrals.
Possible analytic-regulators can be introduced in the propagator power list "FPdenPowers" with the symbol named "apr" by default (whose 
value 0 gives the wanted limit).
*)

ExpansionByRegionGDRI[paraVarSet_List,monomialprf_,polyfactors_List,PFexponents_List,exkinepowerscalingraw_List,po_Integer,OptionsPattern[{UseHomogeneity-> True,ExpansionSpaceDim-> 4,AnalyticRegulator-> apr,KeepEPRarameterAs-> 1,Verbosity-> True}]]:=Module[{svx,sv,lcmN,exkinepowerscaling,Xregpowerraw,Xregpower,Xregpowerscaling,pfsXreg,redpfsXreg,redpfsXregTaylorExpanded,monoprfSP,LOpowerXregion,DRIntegrandXregion},
(***************************A few preliminary checks of the inputs***********************************************************)
(*check the input form of the first four input arguments which all together define the DR-integral*)
If[Union[Map[Head,paraVarSet]]=!={Symbol},Print["The integration variables must all be of the type \"Symbol\"."];Abort[]];
If[Head[monomialprf]===Plus,Print["The 2-th argument saves a purely monomial global pre-factor of parametric-integrand in integration variables (which can be 1)."];Abort[]];
If[(Length[polyfactors]=!=Length[PFexponents])||(Length[polyfactors]<2),Print["The list of \"non-trivial\" polynomial-factors given as the 3-th argument must has the same length as that of the 4-th argument which saves their respective generic exponents, and should have at least two entries (one of the them can be 1)."];Abort[]];
(*check whether the argument that defines the kinematic-limit is specified in the wanted format.*)
With[{svmonolist=Map[Last,exkinepowerscalingraw]},
If[(Union[(svmonolist/.{x-> 1})]=!={1}),
Print["The power-scaling behaviors of all external-kinematics should be specified such that the right-hand sides of each rule
 is purely a monic power of the unique scaling-parameter."];Abort[];]];
(*identify the scaling(expansion)-parameter from the input "exkinepowerscalingraw".*)
svx=With[{svfind=Complement[Variables[exkinepowerscalingraw/.{Rule-> List}],Variables[First[Transpose[exkinepowerscalingraw/.{Rule-> List}]]]]},
If[Length[svfind]===1,First[svfind],Print["The input specifying the power-scaling behavior of external-kinematics is not given in the proper format, and the scaling parameter is not found."];Abort[]]];
If[Not[SubsetQ[{Dr,OptionValue[AnalyticRegulator]},Variables[PFexponents]]],Print["The list of polynomial powers contain symbols other than the specified analytic-regulator (its default name is \"apr\" set in the option 
\"AnalyticRegulator\")!"];Abort[]];
If[OptionValue[Verbosity]===True,Print["The scaling-parameter appearing in the input is : ",svx]];
(**************************Get the region-vectors from ASY*****************************************************)
(*get the list of power-scaling vectors for all scaleful expansion-regions: feeding two products made out of the input-polynomials to WilsonExpand[].*)
Xregpowerraw=WilsonExpand[Apply[Times,Most[polyfactors]],Last[polyfactors],paraVarSet,exkinepowerscalingraw,Delta->OptionValue[UseHomogeneity]];
If[(Xregpowerraw==={})||Not[FreeQ[Xregpowerraw,{}]],Print["No region is found, and the expansion is aborted."];Abort[]];
(* The concrete form of exkinepowerscalingraw itself implies a specific definition of the expansion-parameter, and it can be 
that the region-vectors (power-scaling vector) defined w.r.t this specific expansion-parameter has rational components.
Obtain the integer which is needed to make sure all power-scaling vectors have only integer components, 
which in most considered cases is just 1.*)
lcmN=With[{denfs=DeleteCases[Abs[Flatten[Xregpowerraw]],0]},If[denfs==={},1,Apply[LCM,Map[Denominator,denfs]]]];
If[OptionValue[Verbosity]===True,Print["The least-common-multipler of all denomaintors of region-vectors' components is: ",lcmN]];
(*for the time-being, just rename the expansion-parameter "x" (which is saved in svx) as "sv"*)
exkinepowerscaling=exkinepowerscalingraw/.{Rule[l_,r_]:>Rule[l,r^lcmN] }/.{svx-> sv};
If[OptionValue[Verbosity]===True,Print["The power-scaling of external-kinematics (in dummay x) used now reads as: ",exkinepowerscaling/.{sv-> x}]];
(**************************Apply the expansion-by-regions (Begin)*****************************************************)
Xregpower=If[OptionValue[UseHomogeneity],Map[(#-Min[#])&,(lcmN*Xregpowerraw)],(lcmN*Xregpowerraw)];
If[OptionValue[Verbosity]===True,Print["The list of (reformulated) region-vectors reads as: ",Xregpower]];
Xregpowerscaling=Map[Thread[Rule[paraVarSet,Power[sv,#]]]&,Xregpower];
(*************************************)
pfsXreg=Table[With[{pfTerms=MonomialList[pfor],pfMonomials=(MonomialList[pfor]/.{Times[nc_Integer,rest__]:> Times[1,rest]})},Map[(With[{svPSlist=Map[(Exponent[#,sv,Min])&,pfMonomials/.exkinepowerscaling/.##]},{Min[svPSlist],Total[pfTerms*Power[sv,(svPSlist-Min[svPSlist])]]}])&,Xregpowerscaling]],{pfor,polyfactors}];
(*************************************)
monoprfSP=Map[Exponent[monomialprf,#]&,paraVarSet];
(*"+1" for every integration-variable is to account for the power-scaling contribution of the integration measure*)
LOpowerXregion=Simplify[MapIndexed[((Dot[(monoprfSP+1),#1])+Sum[Times[First[pfsXreg[[pi]][[First[#2]]]],PFexponents[[pi]]],{pi,1,Length[polyfactors]}])&,Xregpower]/.{Dr-> OptionValue[ExpansionSpaceDim],OptionValue[AnalyticRegulator]-> 0}];
If[OptionValue[Verbosity]===True,Print["The LO-power per region: ",LOpowerXregion]];
(*************************************)
(*Now we compose the core piece of the sv-manifested "reduced" form of the integrand (without monomial prefactors of symbolic power-exponent) in each region.*)
redpfsXreg=Table[Product[Power[Last[pfsXreg[[pi]][[n]]],PFexponents[[pi]]],{pi,1,Length[polyfactors]}],{n,1,Length[Xregpower]}];
(*The (po*lcmN-LOpowerXregion[[First[#2]]]) can be negative for some regions, and then Series[] gives 0 as expected.*)
redpfsXregTaylorExpanded=MapIndexed[CompoundExpression[If[OptionValue[Verbosity]===True,Print["The Taylor-expansion-depth w.r.t the LO-result in region-",First[#2]," is ",(po*lcmN-LOpowerXregion[[First[#2]]])+1]],Normal[Series[#1,{sv,0,(po*lcmN-LOpowerXregion[[First[#2]]])}]]]&,redpfsXreg];
(*The following definition agrees with SecDec convention*)
redpfsXregTaylorExpanded=monomialprf*((Power[sv,LOpowerXregion]*redpfsXregTaylorExpanded)//.{sv-> OptionValue[KeepEPRarameterAs]});
(**************************Apply the expansion-by-regions (End)*****************************************************)
Return[{(LOpowerXregion/lcmN),redpfsXregTaylorExpanded}];
];
