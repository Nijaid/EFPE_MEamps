(* ::Package:: *)

BeginPackage["EccentricAmplitudes`Series`",{"EccentricAmplitudes`ConvList`","EccentricAmplitudes`AmplitudeFuncs`"}];


AmplitudesTable::usage = "AmplitudesTable[e,l,m,eval] tabulates the (l,m)-mode amplitudes for an eccentricity e.
If eval is a number, the amplitudes are evaluated numerically at that given precision eval;
if eval is not a number, the amplitudes are returned in symbolic form with eval is the symbol for eccentricity..
{AmplitudesTable, jLimits} are returned, where jLimits is the list of harmonics j}."
recoverSeries::usage = "recoverSeries[e,ell,l,m,eval] returns the harmonic series of the (l,m) mode for a mean anomaly ell.
See AmplitudesTable for the eval option."
recoverConvList::usage = "recoverConvList[e,l,m] returns the ConvList for a specified eccentricity e."


Begin["Private`"]


sumFunc[SeriesFunc_,e_,index_,max_]:=sumFunc[SeriesFunc,e,index,max]=
	Sum[SeriesFunc[e,index,ith],{ith,0,max}];


recoverH2m0[ee_,m_]:=Switch[m,
	0,Return[H200ConvList[[ee]]];,
	2,Return[H220ConvList[[ee]]];,
	_,Message[recoverConvList::nm,m]]
	
recoverConvList::nm="m=`1` is not available";


getConvList[ee_,l_,m_]:=Switch[l,
	2,Return[recoverH2m0[ee,m]];,
	_,Message[recoverConvList::nl,l]]
	
recoverConvList::nl="l=`1` is not available";


recoverConvList[e_,l_,m_]:=Module[{ee=FirstPosition[EccList,_?(e<=#&)]},
	If[MissingQ[ee],Message[recoverConvList::necc,e],
		Return[getConvList[First[ee],l,m]];
	]
]
recoverConvList::necc="e=`1` is not covered by our tables."


N20j0Table[e_,ConvList_,prec_]:=Module[{jList=Range[1,ConvList],tab,jNum,j,jj,jTerm,tabLength},
	jNum=Length[jList];
	tab=ConstantArray[0,2*jNum+1];
	tabLength=Length[tab];
	
	For[jj=1,jj<=jNum,jj++,
		(* Fill in tab from the ends (j=jList[[-1]]) to its center (j=jList[[1]]). *)
		j=jList[[jNum-jj+1]];
		jTerm=N20j0[e,j];
		tab[[jj]]=tab[[tabLength-jj+1]]=If[NumericQ[prec],N[jTerm,prec],jTerm];];
		
	Return[{tab,Range[-ConvList,ConvList]}];]


PFunc=PmW0[#1,2,#2,#3]&;
KFunc=Ks0;

N22j0Table[e_,ConvList_,eval_]:=Module[{jList=Range@@ConvList[[1]],sLimits=ConvList[[2]],jsnmaxs=ConvList[[3]],tab,sList,jj,ss,j,s,jTerm,Pterm,Kterm},
	tab=ConstantArray[0,Length[jList]];
	
	For[jj=1,jj<=Length[jList],jj++,
		j=jList[[jj]];
		sList=Range@@sLimits[[jj]];
		jTerm=0;
		
		For[ss=1,ss<=Length[sList],ss++,
			s=sList[[ss]];
			With[{Pindex=s,Kindex=s-j,Pmax=jsnmaxs[[jj]][[ss]][[1]],Kmax=jsnmaxs[[jj]][[ss]][[2]]},
				Pterm=sumFunc[PFunc,e,Pindex,Pmax];
				Kterm=sumFunc[KFunc,e,Kindex,Kmax];];
			(* Evaluate numerically if prec specifies a precision *)
			jTerm+=If[NumericQ[eval],N[Pterm*Kterm,eval],Pterm*Kterm];
		];
		
		tab[[jj]]=jTerm;
	];
	
	Return[{tab,jList}];
]


AmplitudesTable[ecc_,l_,m_,eval_]:=Module[{ConvList=recoverConvList[ecc,l,m],e},
	(* If prec is a precision, then the amplitudes are evaluated numerically.
	Otherwise, it is used as the symbol for eccentricity. *)
	If[NumericQ[eval],e=ecc,e=eval];
	Switch[m,
		0,N20j0Table[e,ConvList,eval],
		2,N22j0Table[e,ConvList,eval]
	]
]


recoverSeries[e_,ell_,l_,m_,prec_]:=Module[{ampsout=AmplitudesTable[e,l,m,prec],series=0,jj,amps,jList},
	amps=ampsout[[1]];
	jList=ampsout[[2]];
	
	For[jj=1,jj<=Length[jList],jj++,
		series+=amps[[jj]] Exp[-I jList[[jj]] ell];
	];
		
	series
]


End[];


EndPackage[];
