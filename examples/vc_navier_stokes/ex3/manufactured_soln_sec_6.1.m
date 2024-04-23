(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]


(* ::Input::Initialization:: *)

(* Can be used for Dirichlet and traction BCs *)
(* Desired solution *)
u[x_,y_,t_]:=2\[Pi] Cos[2 \[Pi] x]Cos[2 \[Pi] t - 2 \[Pi] y];
v[x_,y_,t_]:= 2\[Pi] Sin[2 \[Pi] x]Sin[2 \[Pi] t - 2 \[Pi] y]-Sin[2 \[Pi] t - 2 \[Pi] x];
p[x_,y_,t_]:=-2\[Pi] Sin[2 \[Pi] t - 2\[Pi] x]Cos[2 \[Pi] t - 2 \[Pi] y];
d[x_,y_]:=0.1-Sqrt[(x-1/2)^2+(y-1/2)^2];
\[Rho][x_,y_,t_]:= \[Rho]0 + (\[Rho]1/2) *(Tanh[d[x, y]/ \[Delta]] + 1);
\[Mu][x_,y_,t_]:= \[Mu]0 + \[Mu]1 +  \[Mu]1 Sin[2 \[Pi] x] Cos[2 \[Pi] y];
(*\[Mu][x_,y_,t_]:= \[Mu]0  + \[Mu]1;*)
(*\[Mu][x_,y_,t_]:= \[Mu]0 + \[Mu]1 + \[Mu]1*Sin[2*\[Pi]*t];*)



(* ::Input::Initialization:: *)
tau11[x_,y_,t_]=  D[u[x,y,t],x];
tau12[x_,y_,t_] = 1/2(D[u[x,y,t],y]+ D[v[x,y,t],x]);
tau22[x_,y_,t_]=  D[v[x,y,t],y];


(* ::Input::Initialization:: *)
(* Divergence free? *)
divu[x_, y_, t_]= (D[u[x,y,t],x]+D[v[x,y,t],y])


(* ::Input::Initialization:: *)
xbulkviscosity[x_,y_,t_] = D[-2/3 \[Mu][x, y, t] divu[x, y, t], x]

ybulkviscosity[x_,y_, t_] = D[-2/3 \[Mu][x, y, t] divu[x, y, t], y]


(* ::Input::Initialization:: *)

mass = Simplify[D[\[Rho][x,y,t], t] + D[\[Rho][x,y,t]u[ x,y,t],x] + D[\[Rho][ x,y,t]v[x,y,t],y]]


(* ::Input::Initialization:: *)
(* Required forcing function *)
fx[x_,y_,t_]= FullSimplify[(D[\[Rho][x,y,t]u[x,y,t],t] + D[\[Rho][x,y,t]u[x,y,t]u[x,y,t],x] +  D[\[Rho][x,y,t]v[x,y,t]u[x,y,t],y]) + D[p[x,y,t],x] - 2(D[\[Mu][x,y,t] tau11[x,y,t],x]+D[\[Mu][x,y,t] tau12[x,y,t],y])-xbulkviscosity[x, y,t]]
fy[x_,y_,t_]= FullSimplify[(D[\[Rho][x,y,t]v[x,y,t],t] + D[\[Rho][x,y,t]u[x,y,t]v[x,y,t],x] +  D[\[Rho][x,y,t]v[x,y,t]v[x,y,t],y]) + D[p[x,y,t],y] - 2(D[\[Mu][x,y,t] tau12[x,y,t],x]+D[\[Mu][x,y,t]tau22[x,y,t],y])-ybulkviscosity[x, y, t]]



