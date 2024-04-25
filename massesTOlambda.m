(* ::Package:: *)

(* Formula to write Sin(beta) and Cos(beta) in terms of Tan(beta) *)
(* Beta chosen as 0<Beta<Pi/2 *)


TOtanb[exp_]:=exp/.{SB->TB/Sqrt[1+TB^2],CB->1/Sqrt[1+TB^2]};


(* ::Input::Initialization:: *)
(* writing \[Lambda]1...\[Lambda]5 in terms of masses, m12^2, mixing angle, beta angle, EW vev, \[Lambda]6, \[Lambda]7*) 


(*  (GUNION and HABER hep-ph/0207010) - Equations D13-D17 *)


(* ::Input::Initialization:: *)
\[Lambda]1[MHH_,Mh0_,M122_,CA_,SA_,SB_,CB_,v_,\[Lambda]6_,\[Lambda]7_]:=(MHH^2*CA^2+Mh0^2*SA^2-M122*SB/CB)/(v^2*CB^2)-3/2*\[Lambda]6*SB/CB+1/2 \[Lambda]7*(SB/CB)^3;


(* ::Input::Initialization:: *)
\[Lambda]2[MHH_,Mh0_,M122_,CA_,SA_,SB_,CB_,v_,\[Lambda]6_,\[Lambda]7_]:= (MHH^2*SA^2+Mh0^2*CA^2-M122*CB/SB)/(v^2*SB^2)+1/2*\[Lambda]6*(CB/SB)^3-3/2 \[Lambda]7*(CB/SB);


(* ::Input::Initialization:: *)
\[Lambda]3[MHH_,Mh0_,MHp_,M122_,CA_,SA_,SB_,CB_,v_,\[Lambda]6_,\[Lambda]7_]:=((MHH^2-Mh0^2)*SA*CA+2*MHp^2*SB*CB-M122)/(v^2*SB*CB)-1/2*\[Lambda]6*CB/SB-1/2 \[Lambda]7*(SB/CB);


\[Lambda]4[MA0_,MHp_,M122_,SB_,CB_,v_,\[Lambda]6_,\[Lambda]7_]:=
((MA0^2-2*MHp^2)SB*CB+M122)/(v^2*SB*CB)-1/2*\[Lambda]6*CB/SB-1/2 \[Lambda]7*(SB/CB);


\[Lambda]5[MA0_,m122_,SB_,CB_,v_,\[Lambda]6_,\[Lambda]7_]:=
(m122-MA0^2*SB*CB)/(v^2*SB*CB)-1/2*\[Lambda]6*CB/SB-1/2 \[Lambda]7*(SB/CB);


(* m122 in terms of lambda5 *)


m122N[MA0_,SB_,CB_,v_,\[Lambda]5_,\[Lambda]6_,\[Lambda]7_]:=
MA0^2*SB*CB+(\[Lambda]5+1/2*\[Lambda]6*CB/SB+1/2 \[Lambda]7*(SB/CB))*(v^2*SB*CB);
