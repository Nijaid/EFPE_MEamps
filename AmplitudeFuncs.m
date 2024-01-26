(* ::Package:: *)

BeginPackage["EccentricAmplitudes`AmplitudeFuncs`"];


N20j0::usage = "N20j0[e,j] gives the 0PN N^{20}_j(e) coefficient."

(*KC::usage = "KC[e,s,k] gives the K^C_{s,k} coefficient"
KS::usage = "KS[e,s,k] gives the K^S_{s,k} coefficient"*)
Ks0::usage = "Ks0[e,s,k] gives the 0PN K^{22}_{s,k}(e) coefficient."
PmW0::usage = "PmW0[e,m,j,k] gives the 0PN P^{m2}_{j,k}(e) coefficient."


Begin["`Private`"]


DBesselJ[n_,z_]:=(BesselJ[n-1,z]-BesselJ[n+1,z])/2;
\[Beta]\[Phi][e\[Phi]_]=(1-Sqrt[1-e\[Phi]^2])/e\[Phi];
vofu0[u_,e_]=u+2ArcTan[Simplify[\[Beta]\[Phi][e] Sin[u]/(1-\[Beta]\[Phi][e] Cos[u])]];


(* K^{22}_{s,k} *)
KC[e_,0,0]:=2Sqrt[1-e^2];
KC[e_,0,k_]:=0;
KC[e_,s_,0]:=2(1-e^2)(4BesselJ[s,s e]+2*If[s>1,Sum[BesselJ[i,i e]BesselJ[s-i,(s-i)e],{i,1,s-1}],0])-2BesselJ[s,s e];
KC[e_,s_,k_]:=KC[e,s,k]=8(1-e^2)BesselJ[s+k,(s+k)e]BesselJ[k,k e];
KS[e_,s_,0]:=4e Sqrt[1-e^2](DBesselJ[s,s e]+If[s>1,Sum[BesselJ[m,m e]DBesselJ[s-m,(s-m)e],{m,1,s-1}],0]);
KS[e_,s_,k_]:=KS[e,s,k]=4e Sqrt[1-e^2](-BesselJ[s+k,(s+k)e]DBesselJ[k,k e]+BesselJ[k,k e]DBesselJ[s+k,(s+k)e]);
Ks0[e_,0,k_]:=KC[e,0,k];
Ks0[e_,s_,k_]:=1/2 (KC[e,Abs[s],k]+Sign[s]KS[e,Abs[s],k]);


(* P^{2W}_{s,k} *)
\[Epsilon]u0[e_,0,0]:=1;
\[Epsilon]u0[e_,k_,s_]:=\[Epsilon]u0[e,k,s]=If[k>=0,If[s==0,-(1/2)e KroneckerDelta[k,1],k/s BesselJ[s-k,s e]], \[Epsilon]u0[e,-k,-s]];
\[CapitalEpsilon]0[e_,j_,k_]:=If[k==0,(-\[Beta]\[Phi][e])^j,
	Block[{ecc,kvar},
	Binomial[kvar-1,kvar-j]Hypergeometric2F1[-j,kvar,kvar-j+1,\[Beta]\[Phi][ecc]^2]\[Beta]\[Phi][ecc]^(kvar-j)/.kvar->k/.ecc->e]];
(* Using Block to avoid evaluations of indeterminate values *)

\[Epsilon]v0[e_,0,0,0]:=1;
\[Epsilon]v0[e_,0,0,ith_]:=0;
\[Epsilon]v0[e_,j_,s_,ith_]:=\[Epsilon]v0[e,j,s,ith]=If[j>=0,\[CapitalEpsilon]0[e,j,ith]\[Epsilon]u0[e,ith,s],
\[Epsilon]v0[e,-j,-s,ith]];

PmW0[e_,m_,j_,k_]:=PmW0[e,m,j,k]=\[Epsilon]v0[e,m,j+m,k]


N20j0[e_,j_]:=Sqrt[2/3]BesselJ[j,j e];


End[];


EndPackage[];
