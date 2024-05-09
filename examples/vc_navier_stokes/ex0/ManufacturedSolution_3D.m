(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]


(* ::Input::Initialization:: *)

(* Can be used for Dirichlet and traction BCs *)
(* Desired solution *)
u[x_,y_,z_,t_]:=2Pi Cos[2 Pi x] Cos[2 Pi y-2 Pi t]Sin[2Pi z];
v[x_,y_,z_,t_]:= 2 Pi Sin[2 Pi x]Sin[2Pi y-2 Pi t]Cos[2Pi z];
w[x_,y_,z_,t_]:=Pi Cos[2 Pi t-2 Pi y] Sin[2 Pi x] (-2 Cos[2 Pi z]-2 Sin[2 Pi z])
p[x_,y_,z_,t_]:=2Pi Sin[2Pi x-2Pi t]Cos[2Pi y-2Pi t]Sin[2Pi z-2Pi t];
d[x_,y_,z_]:=0.1-Sqrt[(x-0.5)^2+(y-0.5)^2 + (z-0.5)^2];
\[Rho][x_,y_,z_,t_]:=\[Rho]0+\[Rho]1/2(1+Tanh[d[x,y,z]/\[Delta]]);
\[Mu][x_,y_,z_,t_]:=\[Mu]1 Sin[2 Pi x]Cos[2 Pi y]Sin[2Pi z]+\[Mu]1 +\[Mu]0;



(* ::Input::Initialization:: *)
tau11[x_,y_,z_]=  D[u[x,y,z,t],x];
tau12[x_,y_,z_] = 1/2(D[u[x,y,z,t],y]+ D[v[x,y,z,t],x]);
tau13[x_,y_,z_] = 1/2(D[u[x,y,z,t],z]+ D[w[x,y,z,t],x]);
tau22[x_,y_,z_]=  D[v[x,y,z,t],y];
tau23[x_,y_,z_] = 1/2(D[v[x,y,z,t],z]+ D[w[x,y,z,t],y]);
tau33[x_,y_,z_]=  D[w[x,y,z,t],z];


(* ::Input::Initialization:: *)
(* Divergence free? *)
D[u[x,y,z,t],x]+D[v[x,y,z,t],y]+D[w[x,y,z,t],z]//FullSimplify


(* ::Input::Initialization:: *)
(* Required forcing function *)
fx[x_,y_,z_,t_]= \[Rho][x,y,z,t](D[u[x,y,z,t],t] + u[x,y,z,t]D[u[x,y,z,t],x] + v[x,y,z,t] D[u[x,y,z,t],y]+w[x,y,z,t] D[u[x,y,z,t],z]) + D[p[x,y,z,t],x] - 2(D[\[Mu][x,y,z,t] tau11[x,y,z],x]+D[\[Mu][x,y,z,t] tau12[x,y,z],y]+D[\[Mu][x,y,z,t] tau13[x,y,z],z])//Simplify
fy[x_,y_,z_,t_]= \[Rho][x,y,z,t](D[v[x,y,z,t],t] + u[x,y,z,t]D[v[x,y,z,t],x] + v[x,y,z,t] D[v[x,y,z,t],y]+w[x,y,z,t] D[v[x,y,z,t],z]) + D[p[x,y,z,t],y] - 2(D[\[Mu][x,y,z,t] tau12[x,y,z],x]+D[\[Mu][x,y,z,t] tau22[x,y,z],y]+D[\[Mu][x,y,z,t]tau23[x,y,z],z])//Simplify
fz[x_,y_,z_,t_]= \[Rho][x,y,z,t](D[w[x,y,z,t],t] + u[x,y,z,t]D[w[x,y,z,t],x] + v[x,y,z,t] D[w[x,y,z,t],y]+w[x,y,z,t] D[w[x,y,z,t],z]) + D[p[x,y,z,t],z] - 2(D[\[Mu][x,y,z,t] tau13[x,y,z],x]+D[\[Mu][x,y,z,t] tau23[x,y,z],y]+D[\[Mu][x,y,z,t] tau33[x,y,z],z])//Simplify


(* ::Input::Initialization:: *)
(*Periodic forces?*)
fx[0,y,z,t]-fx[1,y,z,t]//FullSimplify
fy[0,y,z,t]-fy[1,y,z,t]//FullSimplify
fz[0,y,z,t]-fz[1,y,z,t]//FullSimplify
fx[x,0,z,t]-fx[x,1,z,t]//FullSimplify
fy[x,0,z,t]-fy[x,1,z,t]//FullSimplify
fz[x,0,z,t]-fz[x,1,z,t]//FullSimplify
fx[x,y,0,t]-fx[x,y,1,t]//FullSimplify
fy[x,y,0,t]-fy[x,y,1,t]//FullSimplify
fz[x,y,0,t]-fz[x,y,1,t]//FullSimplify


(* ::Input::Initialization:: *)
Manipulate[ContourPlot[u[x,y,t],{x,0,1},{y,0,1},Contours->10],{t,0,0.1}]


(* ::Input::Initialization:: *)
D[p[x,y,t],x]//FullSimplify


(* ::Input::Initialization:: *)
D[p[x,y,t],y]//FullSimplify


(* ::Input::Initialization:: *)
CellPrint[[HoldForm[fx[x,y,t]],"Output",AutoMultiplicationSymbol->True]]


(* ::Input::Initialization:: *)
fz[x,y,z,t]//InputForm


(* ::Input::Initialization:: *)
\[Rho]0 = 1.0;
\[Rho]1 = 1.0+10;
\[Mu]0 = 0.001;
\[Mu]1 = -\[Mu]0+1;
\[Delta] = 0.01;
Manipulate[ContourPlot[\[Mu][x,y,t],{x,0,1},{y,0,1},Contours->10,PlotLegends->Automatic],{t,0,0.1}]


(* ::Input:: *)
(*(* Tangential tractions at top at bottom y faces*)*)


(* ::Input::Initialization:: *)
\[Mu][x,y,z,t](D[u[x,y,z,t],y]+ D[v[x,y,z,t],x])


(* ::Input::Initialization:: *)
\[Mu][x,y,z,t](D[w[x,y,z,t],y]+ D[v[x,y,z,t],z])//InputForm


(* ::Input:: *)
(*(* Normal tractions at the top and bottom y faces *)*)


(* ::Input::Initialization:: *)
-p[x,y,z,t]+2\[Mu][x,y,z,t]D[v[x,y,z,t],y]//InputForm
