(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]


(* ::Input::Initialization:: *)

(* Can be used for Dirichlet and traction BCs *)
(* Desired solution *)
u[x_,y_,t_]:=2Pi Cos[2 Pi x] Cos[2 Pi y-2 Pi t];
v[x_,y_,t_]:= 2 Pi Sin[2 Pi x]Sin[2Pi y-2 Pi t] + Sin[2 Pi x-2Pi t];
p[x_,y_,t_]:=2Pi Sin[2Pi x-2Pi t]Cos[2Pi y-2Pi t];
d[x_,y_]:=0.1-Sqrt[(x-0.5)^2+(y-0.5)^2];
\[Rho][x_,y_,t_]:=\[Rho]0+\[Rho]1/2(1+Tanh[d[x,y]/\[Delta]]);
\[Mu][x_,y_,t_]:=\[Mu]1 Sin[2 Pi x]Cos[2 Pi y]+\[Mu]1 +\[Mu]0;



(* ::Input::Initialization:: *)
tau11[x_,y_,t_]=  D[u[x,y,t],x];
tau12[x_,y_,t_] = 1/2(D[u[x,y,t],y]+ D[v[x,y,t],x]);
tau22[x_,y_,t_]=  D[v[x,y,t],y];


(* ::Input::Initialization:: *)
(* Divergence free? *)
D[u[x,y,t],x]+D[v[x,y,t],y]


(* ::Input::Initialization:: *)
(* Required forcing function *)
fx[x_,y_,t_]= \[Rho][x,y,t](D[u[x,y,t],t] + u[x,y,t]D[u[x,y,t],x] + v[x,y,t] D[u[x,y,t],y]) + D[p[x,y,t],x] - 2(D[\[Mu][x,y,t] tau11[x,y,t],x]+D[\[Mu][x,y,t] tau12[x,y,t],y])
fy[x_,y_,t_]= \[Rho][x,y,t](D[v[x,y,t],t] + u[x,y,t]D[v[x,y,t],x] + v[x,y,t] D[v[x,y,t],y]) + D[p[x,y,t],y] - 2(D[\[Mu][x,y,t] tau12[x,y,t],x]+D[\[Mu][x,y,t]tau22[x,y,t],y])


(* ::Input::Initialization:: *)
(*Periodic forces?*)
fx[0,y,t]-fx[1,y,t]
fy[0,y,t]-fy[1,y,t]
fx[x,0,t]-fx[x,0,t]
fy[x,1,t]-fy[x,1,t]


(* ::Input::Initialization:: *)
Manipulate[ContourPlot[v[x,y,t],{x,0,1},{y,0,1},Contours->10],{t,0,0.001}]


(* ::Input::Initialization:: *)
\[Rho]0 = 1.0;
\[Rho]1 = -1.0+1000;
\[Delta] = 0.01;
ContourPlot[\[Rho][x,y,0],{x,0,1},{y,0,1},Contours->20]
\[Rho][0.5,0.5,0]/\[Rho][1,1,0]


(* ::Input::Initialization:: *)
D[p[x,y,t],x]//FullSimplify


(* ::Input::Initialization:: *)
D[p[x,y,t],y]//FullSimplify


(* ::Input::Initialization:: *)
CellPrint[[HoldForm[fx[x,y,t]],"Output",AutoMultiplicationSymbol->True]]


(* ::Input::Initialization:: *)
fy[x,y,t]//InputForm


(* ::Input::Initialization:: *)
\[Rho]0 = 1.0;
\[Rho]1 = 1.0+10;
\[Mu]0 = 0.001;
\[Mu]1 = -\[Mu]0+1;
\[Delta] = 0.01;
Manipulate[ContourPlot[fy[x,y,t],{x,0,1},{y,0,1},Contours->10,PlotLegends->Automatic],{t,0,0.1}]


(* ::Input::Initialization:: *)
(*Tangential traction at top and bottom face*)
\[Mu][x,y,t](D[u[x,y,t],y]+ D[v[x,y,t],x])


(* ::Input::Initialization:: *)
(*Normal traction at top and bottom face*)
-p[x,y,t]+2\[Mu][x,y,t]D[v[x,y,t],y]//InputForm


(* ::Input::Initialization:: *)
trac[x_,y_,t_]=-p[x,y,t]+2\[Mu][x,y,t]D[v[x,y,t],y]


(* ::Input::Initialization:: *)
trac[0,1,t]-trac[1,1,t]


(* ::Input:: *)
(**)
